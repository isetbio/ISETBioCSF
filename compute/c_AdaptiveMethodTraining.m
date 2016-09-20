function c_EffectOfTrainingSize(varargin)
% c_AdaptiveMethodTraining(varargin)
%
% Explore learning classifier in the context of adaptive psychophysical
% methods.
%
% This looks at L+M detection thrsholds, by default at a moderate spatial frequency (10
% cpd).  The stimulus size is inversely proportional to spatial frequency,
% so using a larger spatial frequency speeds things up.
%
% Key/value pairs
%   'nResponseSamples' - value (default 1000).  Number of precomputed responses to
%       draw from.
%   'trainingFraction' - value (0.8). Fraction of responses to use for
%       training.  Rest are used for test
%   'nAdaptiveTrialsPerStaircase' - value (default 400). Number of trials to run the adaptive
%       procedure
%   'nTestTrialsPerContrast' - value (default 20).  Number of test trials per contrast,
%       after training.
%   'cyclesPerDegree' - value (default 10). Spatial frequency of grating to be investigated.
%       This is reciprocally related to the size of the image/mosaic.
%   'computeResponses' - true/false (default false).  Compute responses.
%   'plotPsychometric' - true/false (default true).  Plot psychometric
%       functions.
%   'plotTrainingSize' - true/false (default true).  Plot results.

%% Parse input
p = inputParser;
p.addParameter('nResponseSamples',1000,@isnumeric);
p.addParameter('trainingFraction',0.8,@isnumeric);
p.addParameter('nAdaptiveTrialsPerStaircase',400,@isnumeric);
p.addParameter('nTestTrialsPerContrast',20,@isnumeric);
p.addParameter('cyclesPerDegree',10,@isnumeric);
p.addParameter('computeResponses',false,@islogical);
p.addParameter('fitPsychometric',true,@islogical);
p.addParameter('plotPsychometric',true,@islogical);
p.addParameter('plotTrainingSize',true,@islogical);
p.parse(varargin{:});

%% Get the parameters we need
%
% Start with default
rParams = colorGaborResponseParamsGenerate;

% Get stimulus parameters correct
%
% The stimulus was half-cosine windowed to contain 7.5 cycles.  We set
% our half-cosine window to match that and also make the field of view
% just a tad bigger.
rParams.gaborParams.windowType = 'halfcos';
rParams.gaborParams.cyclesPerDegree = p.Results.cyclesPerDegree;
rParams.gaborParams.gaussianFWHMDegs = 3.75*(1/rParams.gaborParams.cyclesPerDegree);
rParams.gaborParams.fieldOfViewDegs = 2.1*rParams.gaborParams.gaussianFWHMDegs;

% Keep mosaic size in lock step with stimulus.  This is also forced before
% the mosaic is created, but we need it here so that filenames are
% consistent.  It is possible that we should not have a separate mosaic
% size field, and just alwasy force it to match the scene.
rParams.mosaicParams.fieldOfViewDegs = rParams.gaborParams.fieldOfViewDegs;

% Set background luminance
%
% We start with a base luminance that we know is about mid-gray on the
% monitor we specify.  To change luminance, we specify a scale factor.
% This is eventually applied both to the background luminance and to the
% monitor channel spectra, so that we don't get unintersting out of gamut errors.
baseLum = 50;
theLum = 340;
rParams.backgroundParams.backgroundxyY = [0.33 0.33 baseLum]';
rParams.backgroundParams.monitorFile = 'CRT-MODEL';
rParams.backgroundParams.leakageLum = 1.0;
rParams.backgroundParams.lumFactor = theLum/baseLum;

% Pupil size.  They used a 2mm artificial pupil
oiParams.pupilDiamMm = 2;

% Set duration equal to sampling interval to do just one frame.
%
% Their intervals were 100 msec each.
rParams.temporalParams.simulationTimeStepSecs = 100/1000;
rParams.temporalParams.stimulusDurationInSeconds = rParams.temporalParams.simulationTimeStepSecs;
rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.simulationTimeStepSecs;
rParams.temporalParams.secondsToInclude = rParams.temporalParams.simulationTimeStepSecs;

% Their main calculation was without eye movements
rParams.temporalParams.eyesDoNotMove = true;

% Set up mosaic parameters for just one stimulus time step
rParams.mosaicParams.timeStepInSeconds = rParams.temporalParams.simulationTimeStepSecs;
rParams.mosaicParams.integrationTimeInSeconds = rParams.mosaicParams.timeStepInSeconds;
rParams.mosaicParams.isomerizationNoise = true;
rParams.mosaicParams.osNoise = true;
rParams.mosaicParams.osModel = 'Linear';

%% Parameters that define the LM instances we'll generate here
%
% Use default LMPlane.
testDirectionParams = LMPlaneInstanceParamsGenerate;
testDirectionParams.startAngle = 45;
testDirectionParams.deltaAngle = 90;
testDirectionParams.nAngles = 1;

% Number of contrasts to run in each color direction
testDirectionParams.nContrastsPerDirection = 20;
testDirectionParams.lowContrast = 0.0001;
testDirectionParams.highContrast = 0.1;
testDirectionParams.contrastScale = 'log';

%% Parameters related to how we find thresholds from responses
%
% Use default
thresholdParams = thresholdParamsGenerate;

%% Loop over trianing samples

% Set number of trials
testDirectionParams.trialsNum = p.Results.nResponseSamples;
nTrainingSamples = round(p.Results.trainingFraction*p.Results.nResponseSamples);

%% Compute response instances, which serve as training samples
if (p.Results.computeResponses)
    t_colorGaborConeCurrentEyeMovementsResponseInstances('rParams',rParams,'testDirectionParams',testDirectionParams,'compute',true,'visualizeResponses',false);
end

%% Read/write object
rwObject = IBIOColorDetectReadWriteBasic;
readProgram = 't_colorGaborConeCurrentEyeMovementsResponseInstances';
writeProgram = mfilename;

%% Read training samples for the no stimulus condition
fprintf('Reading no stimulus data ... ');
colorModulationParamsTemp = rParams.colorModulationParams;
colorModulationParamsTemp.coneContrasts = [0 0 0]';
colorModulationParamsTemp.contrast = 0;
paramsList = {rParams.gaborParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
noStimData = rwObject.read('responseInstances',paramsList,readProgram);
ancillaryData = rwObject.read('ancillaryData',paramsList,readProgram);
fprintf('done\n');
responseDimension = length(noStimData.responseInstanceArray(1).theMosaicIsomerizations(:));
testConeContrasts = ancillaryData.testConeContrasts;
if (size(testConeContrasts,2) ~= 1)
    error('We assume only one vector of testConeContrasts, but there is more than one');
end
testContrasts = ancillaryData.testContrasts;

%% Read in set of precomputed training samples for each contrast
stimData = cell(testDirectionParams.nContrastsPerDirection,1);
for cc = 1:testDirectionParams.nContrastsPerDirection
    fprintf('Reading stimulus data for contrast %d of %d ...',cc,testDirectionParams.nContrastsPerDirection);
    colorModulationParamsTemp = rParams.colorModulationParams;
    colorModulationParamsTemp.coneContrasts = testConeContrasts;
    colorModulationParamsTemp.contrast = testContrasts(cc);
    paramsList = {rParams.gaborParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
    stimData{cc} = rwObject.read('responseInstances',paramsList,readProgram);
    fprintf('done\n');
    if (numel(stimData{cc}.responseInstanceArray) ~= testDirectionParams.trialsNum)
        error('Inconsistent number of trials');
    end
end

%% Set up staircases
%
% Some parameters
numTrialsPerStaircase = p.Results.nAdaptiveTrialsPerStaircase;
maxContrast = testDirectionParams.highContrast;
minContrast = testDirectionParams.lowContrast;
stepSizes = [0.05 0.01 0.005 0.001 0.0005 0.0001];
nUps = [3 2 3];
nDowns = [1 1 2];
nInterleavedStaircases = length(nUps);

% Initialize staircase objects
for k = 1:nInterleavedStaircases;
    st{k} = Staircase('standard',maxContrast, ...
        'StepSizes', stepSizes, 'NUp', nUps(k), 'NDown', nDowns(k), ...
        'MaxValue', maxContrast, 'MinValue', minContrast);
end

%% Simulate interleaved staircases and learn decision rule
%
% Make space to hold sensor responses
nTotalTrials = nInterleavedStaircases*numTrialsPerStaircase;
nTrials0 = 0; nTrials1 = 0;
trialResponseData0 = zeros(2*responseDimension,nTotalTrials);
trialResponseData1	 = zeros(2*responseDimension,nTotalTrials);
for i = 1:numTrialsPerStaircase
    staircaseOrder = randperm(nInterleavedStaircases);
    for k = 1:nInterleavedStaircases
        
        % Find contrast to run, making sure it is one we have responses for
        comparisonContrast = getCurrentValue(st{staircaseOrder(k)});
        if (comparisonContrast > maxContrast)
            comparisonContrast = maxContrast;
        elseif (comparisonContrast < minContrast)
            comparisonContrast = minContrast;
        else
            contrastDiffs = abs(testContrasts-comparisonContrast);
            [~,index] = sort(contrastDiffs);
            if (index(1) > 1)
                if (round(rand))
                    index(1) = index(1) - 1;
                end
            end
            comparisonContrast = testContrasts(index(1));
        end
        contrastIndex = find(testContrasts == comparisonContrast);
        if (isempty(contrastIndex))
            error('Oops.  Comparsion contrast not in test set');
        end
        
        % Simulate trial
        noiseIntervalResponseIndex = randi(nTrainingSamples);
        noiseIntervalResponse = noStimData.responseInstanceArray(noiseIntervalResponseIndex).theMosaicIsomerizations(:);
        signalIntervalResponseIndex = randi(nTrainingSamples);
        signalIntervalResponse = stimData{contrastIndex}.responseInstanceArray(signalIntervalResponseIndex).theMosaicIsomerizations(:);
        
        % Determine which trial type (noise-signal, signal-noise) this is,
        % and construct trial response
        trialType = round(rand);
        if (trialType == 0)
            trialResponse = [noiseIntervalResponse ; signalIntervalResponse];
        else
            trialResponse = [signalIntervalResponse ; noiseIntervalResponse];
        end
        
        % Classify trial response.  Just guess if we haven't seen exemplars
        % of each type.  Use nearest neighbor classification on trials we
        % have seen if we have.
        if (nTrials0 == 0 & nTrials1 == 0)
            responseType = round(rand);
        else
            trialDistances0 = zeros(nTrials0,1);
            for ii = 1:nTrials0
                trialDistances0(ii) = norm(trialResponse-trialResponseData0(:,ii));
            end
            trialDistances1 = zeros(nTrials1,1);
            for ii = 1:nTrials1
                trialDistances1(ii) = norm(trialResponse-trialResponseData1(:,ii));
            end
            if (min(trialDistances0) < min(trialDistances1))
                responseType = 0;
            else
                responseType = 1;
            end
        end
        
        % Decide whether response was correct or not
        if (responseType == trialType)
            response = 1;
        else
            response = 0;
        end
        
        % Save trial response data
        if (trialType == 0)
            nTrials0 = nTrials0+1;
            trialResponseData0(:,nTrials0) = trialResponse;
        else
            nTrials1 = nTrials1+1;
            trialResponseData1(:,nTrials1) = trialResponse;
        end
        
        % Update staircase
        st{staircaseOrder(k)} = updateForTrial(st{staircaseOrder(k)},comparisonContrast,response);
        
        % Report
        fprintf('Simulating staircase %d, trial %d, contrast %0.5f, trial type = %d, response type = %d, response %d\n',staircaseOrder(k),i,comparisonContrast,trialType,responseType,response);      
    end
end

%% Now freeze decision rule and measure psychometric function
nCorrects = zeros(length(testContrasts),1);
for cc = 1:length(testContrasts)    
    for jj = 1:p.Results.nTestTrialsPerContrast; 
        
        % Simulate trial.  These trials come from test set of sampled
        % responses.ß
        noiseIntervalResponseIndex = randi(p.Results.nResponseSamples-nTrainingSamples);
        noiseIntervalResponse = noStimData.responseInstanceArray(nTrainingSamples+noiseIntervalResponseIndex).theMosaicIsomerizations(:);
        signalIntervalResponseIndex = randi(p.Results.nResponseSamples-nTrainingSamples);
        signalIntervalResponse = stimData{cc}.responseInstanceArray(nTrainingSamples+signalIntervalResponseIndex).theMosaicIsomerizations(:);
        
        % Determine which trial type (noise-signal, signal-noise) this is,
        % and construct trial response
        trialType = round(rand);
        if (trialType == 0)
            trialResponse = [noiseIntervalResponse ; signalIntervalResponse];
        else
            trialResponse = [signalIntervalResponse ; noiseIntervalResponse];
        end
        
        % Classify trial response.  Just guess if we haven't seen exemplars
        % of each type.  Use nearest neighbor classification on trials we
        % have seen if we have.
        if (nTrials0 == 0 & nTrials1 == 0)
            responseType = round(rand);
        else
            trialDistances0 = zeros(nTrials0,1);
            for ii = 1:nTrials0
                trialDistances0(ii) = norm(trialResponse-trialResponseData0(:,ii));
            end
            trialDistances1 = zeros(nTrials1,1);
            for ii = 1:nTrials1
                trialDistances1(ii) = norm(trialResponse-trialResponseData1(:,ii));
            end
            if (min(trialDistances0) < min(trialDistances1))
                responseType = 0;
            else
                responseType = 1;
            end
        end
        
        % Decide whether response was correct and keep track
        if (responseType == trialType)
            nCorrects(cc) = nCorrects(cc)+1;
        end
    end
end

%% Analyze staircase data
valuesStair = []; responsesStair = [];
nTrialsKeep = 100;
nTrialsDiscard = numTrialsPerStaircase-nTrialsKeep;
for k = 1:nInterleavedStaircases
    threshStair(k) = getThresholdEstimate(st{k});
    [valuesSingleStair{k},responsesSingleStair{k}] = getTrials(st{k},nTrialsDiscard);
    valuesStair = [valuesStair valuesSingleStair{k}];
    responsesStair = [responsesStair responsesSingleStair{k}];
end
[meanValues,nCorrectStair,nTrialsStair] = GetAggregatedStairTrials(valuesStair,responsesStair,20);

%% Fit staircase data using Palamedes
paramsFree     = [1 1 0 0];
PF             = @PAL_Weibull;

% Some optimization settings for the fit
options             = optimset('fminsearch');
options.TolFun      = 1e-09;
options.MaxFunEvals = 1000;
options.MaxIter     = 1000;
options.Display     = 'off';

% Search grid specification for Palemedes
gridLevels = 100;
searchGrid.alpha = logspace(log10(testContrasts(1)),log10(testContrasts(end)),gridLevels);
searchGrid.beta = 10.^linspace(-2,2,gridLevels);
searchGrid.gamma = 0.5;
searchGrid.lambda = 0.0;

% Use Palamedes grid search method
[paramsValues,LL,flag] = PAL_PFML_Fit(valuesStair(:), responsesStair(:), ones(size(responsesStair(:))), ...
    searchGrid, paramsFree, PF, 'SearchOptions', options);

% Get threshold and deal with catastrophic cases
threshold = PF(paramsValues, thresholdParams.criterionFraction, 'inverse');
if (threshold < 0 | threshold > 1 | ~isreal(threshold) | isinf(threshold))
    threshold = NaN;
end

%% Provide fit psychometric function on passed stimulus levels
% if (~isnan(threshold))
%     fitFractionCorrect = PF(paramsValues,fitStimLevels);
% else
%     fitFractionCorrect = NaN*ones(size(fitStimLevels));
% end

%% Make a figure of what happened during the staircase
stairFig = figure; clf;
colors = ['r' 'g' 'b' 'y' 'c'];
subplot(1,2,1); hold on
for k = 1:nInterleavedStaircases
    xvalues = 1:nTrialsKeep;
    index = find(responsesSingleStair{k} == 0);
    plot(xvalues,log10(valuesSingleStair{k}),[colors(k) '-']);
    plot(xvalues,log10(valuesSingleStair{k}),[colors(k) 'o'],'MarkerFaceColor',colors(k),'MarkerSize',6);
    if (~isempty(index))
        plot(xvalues(index),log10(valuesSingleStair{k}(index)),[colors(k) 'o'],'MarkerFaceColor','w','MarkerSize',6);
    end
end
plot(xvalues,log10(threshold)*ones(1,numTrialsPerStaircase),'k');
xlabel('Trial Number','FontSize',16);
ylabel('Log10 Contrast','FontSize',16);
title(sprintf('TAFC staircase plot'),'FontSize',16);

subplot(1,2,2); hold on
plot(log10(meanValues),nCorrectStair./nTrialsStair,'ko','MarkerSize',6,'MarkerFaceColor','k');
plot(log10([threshold threshold]),[0 thresholdParams.criterionFraction],'r','LineWidth',2);
xlabel('Log10 Contrast','FontSize',16);
ylabel('Prob Correct','FontSize',16);
title(sprintf('TAFC staircase psychometric function'),'FontSize',16);
%xlim([comparisonStimuli(1)-testStimulus comparisonStimuli(end)-testStimulus])
ylim([0 1]);
% if (exist('FigureSave','file'))
%     FigureSave('StaircaseFC',gcf','pdf');
% else
%     saveas(gcf','StaircaseFC','pdf');
% end

%% Fit psychometric functions
if (p.Results.fitPsychometric)
    if (p.Results.computeMLPTThresholds)
        thresholdParams.method = 'mlpt';
        effectOfTrainingSize.mlptThresholds(tt) = t_plotGaborDetectThresholdsOnLMPlane('rParams',rParams,'LMPlaneInstanceParams',testDirectionParams,'thresholdParams',thresholdParams, ...
            'plotPsychometric',p.Results.plotPsychometric,'plotEllipse',false);
        close all;
    end
end


%% Write out the data
%
% Set trialsNum to 0 to define a summary directory name
% fprintf('Writing performance data ... ');
% testDirectionParams.trialsNum = 0;
% paramsList = {rParams.gaborParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
% rwObject = IBIOColorDetectReadWriteBasic;
% writeProgram = mfilename;
% rwObject.write('effectOfTrainingSize',effectOfTrainingSize,paramsList,writeProgram);
% fprintf('done\n');
%
% %% Make a plot of estimated threshold versus training set size
% %
% % The way the plot is coded counts on the test contrasts never changing
% % across the conditions, which we could explicitly check for here.
% if (p.Results.plotTrainingSize)
%     fprintf('Reading performance data ...');
%     testDirectionParams.trialsNum = 0;
%     paramsList = {rParams.gaborParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
%     rwObject = IBIOColorDetectReadWriteBasic;
%     writeProgram = mfilename;
%     effectOfTrainingSize = rwObject.read('effectOfTrainingSize',paramsList,writeProgram);
%     fprintf('done\n');
%
%     hFig = figure; clf; hold on
%     fontBump = 4;
%     set(gca,'FontSize', rParams.plotParams.axisFontSize+fontBump);
%     plot(log10(effectOfTrainingSize.nTrainingSamplesList),[effectOfTrainingSize.mlptThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1), ...
%         'ro','MarkerSize',rParams.plotParams.markerSize,'MarkerFaceColor','r');
%     plot(log10(effectOfTrainingSize.nTrainingSamplesList),[effectOfTrainingSize.mlpeThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1), ...
%         'go','MarkerSize',rParams.plotParams.markerSize,'MarkerFaceColor','g');
%     plot(log10(effectOfTrainingSize.nTrainingSamplesList),[effectOfTrainingSize.svmThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1), ...
%         'bo','MarkerSize',rParams.plotParams.markerSize,'MarkerFaceColor','b');
%     plot(log10(effectOfTrainingSize.nTrainingSamplesList),[effectOfTrainingSize.mlptThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1), ...
%         'r','LineWidth',rParams.plotParams.lineWidth);
%     plot(log10(effectOfTrainingSize.nTrainingSamplesList),[effectOfTrainingSize.mlpeThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1), ...
%         'g','LineWidth',rParams.plotParams.lineWidth);
%     plot(log10(effectOfTrainingSize.nTrainingSamplesList),[effectOfTrainingSize.svmThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1), ...
%         'b','LineWidth',rParams.plotParams.lineWidth);
%     xlabel('Log10 Training Set Size', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
%     ylabel('Threshold Contrast', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
%     xlim([1 4]); ylim([0 1e-2]);
%     legend({'MaxLikelihood, signal known', 'MaxLikelihood, signal estimated', sprintf('SVM, PCA %d',thresholdParams.PCAComponents)},'Location','NorthEast','FontSize',rParams.plotParams.labelFontSize+fontBump);
%     box off; grid on
%     title(sprintf('Effect of Training Set Size, %0.2f deg mosaic',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
%     rwObject.write('effectOfTrainingSize',hFig,paramsList,writeProgram,'Type','figure');
% end
