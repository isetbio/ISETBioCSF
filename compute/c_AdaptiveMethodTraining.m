function c_AdaptiveMethodTraining(varargin)
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
%   'learningMethod' - string (default 'nnmse').  Learning rule to use.
%       'nnmes' - nearest neighbor, squared error metric
%       'nncorr' - nearest neighbor, correlation-based error metric
%       'uniflikely' - compute likelyhood that each interval is uniform field.
%   'filterMethod' - string (default 'none').  Filter the data and if so how?
%       'none' - don't do any filtering
%       'lowpass' - lowpass filtering and downsample.
%   'lowpassParam' - value (default XX).  Parameter to use for lowpass filter
%   'nResponseSamples' - value (default 1000).  Number of precomputed responses to
%       draw from.
%   'trainingFraction' - value (0.8). Fraction of responses to use for
%       training.  Rest are used for test
%   'nTrialsPerStaircase' - value (default 400). Number of trials to run the adaptive
%       procedure
%   'nTrialsPerTestContrast' - value (default 20).  Number of test trials per contrast,
%       after training.
%   'cyclesPerDegree' - value (default 10). Spatial frequency of grating to be investigated.
%       This is reciprocally related to the size of the image/mosaic.
%   'computeResponses' - true/false (default false).  Compute responses.
%   'plotPsychometric' - true/false (default true).  Plot psychometric
%       functions.
%   'plotTrainingSize' - true/false (default true).  Plot results.

%% Parse input
p = inputParser;
p.addParameter('learningMethod','nnmse',@ischar);
p.addParameter('filterMethod','none',@ischar);
p.addParameter('lowpassParam',10,@isnumeric);
p.addParameter('nResponseSamples',1000,@isnumeric);
p.addParameter('trainingFraction',0.6,@isnumeric);
p.addParameter('nTrialsPerStaircase',150,@isnumeric);
p.addParameter('nTrialsPerTestContrast',20,@isnumeric);
p.addParameter('cyclesPerDegree',10,@isnumeric);
p.addParameter('computeResponses',false,@islogical);
p.addParameter('fitPsychometric',true,@islogical);
p.addParameter('plotPsychometric',true,@islogical);
p.addParameter('plotTrainingSize',true,@islogical);
p.parse(varargin{:});

%% Get the parameters we need
%
% Start with default
rParams = responseParamsGenerate;

% Get stimulus parameters correct
%
% The stimulus was half-cosine windowed to contain 7.5 cycles.  We set
% our half-cosine window to match that and also make the field of view
% just a tad bigger.
rParams.spatialParams.windowType = 'halfcos';
rParams.spatialParams.cyclesPerDegree = p.Results.cyclesPerDegree;
rParams.spatialParams.gaussianFWHMDegs = 3.75*(1/rParams.spatialParams.cyclesPerDegree);
rParams.spatialParams.fieldOfViewDegs = 2.1*rParams.spatialParams.gaussianFWHMDegs;

% Keep mosaic size in lock step with stimulus.  This is also forced before
% the mosaic is created, but we need it here so that filenames are
% consistent.  It is possible that we should not have a separate mosaic
% size field, and just alwasy force it to match the scene.
rParams.mosaicParams.fieldOfViewDegs = rParams.spatialParams.fieldOfViewDegs;

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
rParams.temporalParams.stimulusDurationInSeconds = 100/1000;
rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.stimulusDurationInSeconds;
rParams.temporalParams.secondsToInclude = rParams.temporalParams.stimulusDurationInSeconds;

% Their main calculation was without eye movements
rParams.temporalParams.emPathType = 'none';

% Set up mosaic parameters for just one stimulus time step
rParams.mosaicParams.integrationTimeInSeconds = rParams.temporalParams.stimulusDurationInSeconds;
rParams.mosaicParams.isomerizationNoise = true;
rParams.mosaicParams.osNoise = true;
rParams.mosaicParams.osModel = 'Linear';

%% Parameters that define the LM instances we'll generate here
%
% Use default LMPlane.
testDirectionParams = instanceParamsGenerate;
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
    t_coneCurrentEyeMovementsResponseInstances('rParams',rParams,'testDirectionParams',testDirectionParams,'compute',true,'visualizeResponses',false);
end

%% Read/write object
rwObject = IBIOColorDetectReadWriteBasic;
readProgram = 't_coneCurrentEyeMovementsResponseInstances';
writeProgram = mfilename;

%% Read training samples for the no stimulus condition
fprintf('Reading no stimulus data ... ');
colorModulationParamsTemp = rParams.colorModulationParams;
colorModulationParamsTemp.coneContrasts = [0 0 0]';
colorModulationParamsTemp.contrast = 0;
paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
noStimData = rwObject.read('responseInstances',paramsList,readProgram);
ancillaryData = rwObject.read('ancillaryData',paramsList,readProgram);
fprintf(' done\n');
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
    paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
    stimData{cc} = rwObject.read('responseInstances',paramsList,readProgram);
    fprintf(' done\n');
    if (numel(stimData{cc}.responseInstanceArray) ~= testDirectionParams.trialsNum)
        error('Inconsistent number of trials');
    end
end

%% Start setting up learning structure
%
% Curiously, the syntax ancillaryData.theMosaic.pattern causes an
% error, while theMosaic.pattern is fine.  But we want to keep
% theMosaic in learningStructure so it is easy to pass around.
learningStructure.theMosaic = ancillaryData.theMosaic;
theMosaic = learningStructure.theMosaic;
for ii = 1:3
    learningStructure.coneIndex{ii} = find(theMosaic.pattern(:) == ii+1);
end
learningStructure.lowpassParam = p.Results.lowpassParam;
clear theMosaic

%% Let's take a look at the mosaic responses at high contrast
figure; clf;
aNoiseIntervalResponse = noStimData.responseInstanceArray(1).theMosaicIsomerizations;
aSignalIntervalResponse = stimData{end}.responseInstanceArray(1).theMosaicIsomerizations;
aNoiseIntervalResponseScale = zeros(size(aNoiseIntervalResponse));
aSignalIntervalResponseScale = zeros(size(aSignalIntervalResponse));
for ii = 1:3
    aNoiseIntervalResponseScale(learningStructure.coneIndex{ii}) = aNoiseIntervalResponse(learningStructure.coneIndex{ii})/mean(aNoiseIntervalResponse(learningStructure.coneIndex{ii}));
    aSignalIntervalResponseScale(learningStructure.coneIndex{ii}) = aSignalIntervalResponse(learningStructure.coneIndex{ii})/mean(aNoiseIntervalResponse(learningStructure.coneIndex{ii}));
end
[imRows,imCols] = size(aNoiseIntervalResponse);
pixelPerm = randperm(imRows*imCols);
aNoiseIntervalResponseScramble = reshape(aNoiseIntervalResponseScale(pixelPerm),imRows,imCols);
aSignalIntervalResponseScramble = reshape(aSignalIntervalResponseScale(pixelPerm),imRows,imCols);
subplot(2,2,1);
imshow(0.8*aNoiseIntervalResponseScale/mean(aNoiseIntervalResponseScale(:)));
title('Blank');
subplot(2,2,2);
imshow(0.8*aSignalIntervalResponseScale/mean(aNoiseIntervalResponseScale(:)));
title('Grating');
subplot(2,2,3);
imshow(0.8*aNoiseIntervalResponseScramble/mean(aNoiseIntervalResponseScale(:)));
title('Scrambled Blank');
subplot(2,2,4);
imshow(0.8*aSignalIntervalResponseScramble/mean(aNoiseIntervalResponseScale(:)));
title('Scrambled Grating');

%% Set up staircases
%
% Some parameters
numTrialsPerStaircase = p.Results.nTrialsPerStaircase;
maxContrast = testDirectionParams.highContrast;
minContrast = testDirectionParams.lowContrast;
stepSizes = [0.05 0.01];
nUps = [3 2 3];
nDowns = [1 1 2];
nInterleavedStaircases = length(nUps);
if (nInterleavedStaircases*p.Results.nTrialsPerStaircase > nTrainingSamples)
    error('Can only use each training sample once during learning, need more training samples');
end

% Initialize staircase objects
for k = 1:nInterleavedStaircases;
    st{k} = Staircase('standard',maxContrast, ...
        'StepSizes', stepSizes, 'NUp', nUps(k), 'NDown', nDowns(k), ...
        'MaxValue', maxContrast, 'MinValue', minContrast);
end

%% Simulate interleaved staircases and learn decision rule on a trial by trial basis
%
% Make space to hold sensor responses
nTotalTrials = nInterleavedStaircases*numTrialsPerStaircase;
learningStructure.nTrials0 = 0; learningStructure.nTrials1 = 0;
learningStructure.trialResponseData0 = zeros(2*responseDimension,nTotalTrials);
learningStructure.trialResponseData1 = zeros(2*responseDimension,nTotalTrials);
sampleCount = 1;
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
        
        % Simulate trial. Every response sample only gets used only once.
        %
        % First get noise and signal interval responses
        noiseIntervalResponse = noStimData.responseInstanceArray(sampleCount).theMosaicIsomerizations;
        signalIntervalResponse = stimData{contrastIndex}.responseInstanceArray(sampleCount).theMosaicIsomerizations;
        sampleCount = sampleCount+1;
        
        % Filter if desired
        noiseIntervalResponse = filterAndVectorizeResponse(p.Results.filterMethod,noiseIntervalResponse,learningStructure);
        signalIntervalResponse = filterAndVectorizeResponse(p.Results.filterMethod,signalIntervalResponse,learningStructure);
        
        % Then simulate trial
        [response,learningStructure] = simulateTrial(p.Results.learningMethod,true,noiseIntervalResponse,signalIntervalResponse,learningStructure);
        
        % Update staircase
        st{staircaseOrder(k)} = updateForTrial(st{staircaseOrder(k)},comparisonContrast,response);
        
        % Report
        if (rem(i,40) == 0)
            fprintf('Simulating staircase %d, trial %d, contrast %0.5f, response %d\n',staircaseOrder(k),i,comparisonContrast,response);
        end
    end
end

%% Summarize staircase data
valuesStair = []; responsesStair = [];
for k = 1:nInterleavedStaircases
    threshStair(k) = getThresholdEstimate(st{k});
    [valuesSingleStair{k},responsesSingleStair{k}] = getTrials(st{k});
    valuesStair = [valuesStair valuesSingleStair{k}];
    responsesStair = [responsesStair responsesSingleStair{k}];
end
[meanValues,nCorrectStair,nTrialsStair] = GetAggregatedStairTrials(valuesStair,responsesStair,30);

%% Make a figure of what happened during the staircases
stairFig = figure; clf;
colors = ['r' 'g' 'b' 'y' 'c'];
subplot(1,2,1); hold on
for k = 1:nInterleavedStaircases
    xvalues = 1:p.Results.nTrialsPerStaircase;
    index = find(responsesSingleStair{k} == 0);
    plot(xvalues,log10(valuesSingleStair{k}),[colors(k) '-']);
    plot(xvalues,log10(valuesSingleStair{k}),[colors(k) 'o'],'MarkerFaceColor',colors(k),'MarkerSize',6);
    if (~isempty(index))
        plot(xvalues(index),log10(valuesSingleStair{k}(index)),[colors(k) 'o'],'MarkerFaceColor','w','MarkerSize',6);
    end
end
xlabel('Trial Number','FontSize',16);
ylabel('Log10 Contrast','FontSize',16);
title(sprintf('TAFC staircase plot'),'FontSize',16);

subplot(1,2,2); hold on
plot(log10(meanValues),nCorrectStair./nTrialsStair,'ko','MarkerSize',6,'MarkerFaceColor','k');
xlabel('Log10 Contrast','FontSize',16);
ylabel('Prob Correct','FontSize',16);
title(sprintf('TAFC staircase psychometric function'),'FontSize',16);
ylim([0 1]);

%% Now freeze decision rule and simulate out psychometric function
nCorrects = zeros(length(testContrasts),1);
if ((p.Results.nResponseSamples-nTrainingSamples) < length(testContrasts)*p.Results.nTrialsPerTestContrast)
    error('Not enough test samples, need more samples');
end
sampleCount = 1;
for cc = 1:length(testContrasts)
    fprintf('Psychometric function calculation, conrast %d of %d\n',cc,length(testContrasts));
    for jj = 1:p.Results.nTrialsPerTestContrast;
        
        % Simulate trial. We only use each sample once
        noiseIntervalResponse = noStimData.responseInstanceArray(nTrainingSamples+sampleCount).theMosaicIsomerizations;
        signalIntervalResponse = stimData{cc}.responseInstanceArray(nTrainingSamples+sampleCount).theMosaicIsomerizations;
        sampleCount = sampleCount+1;
        
        % Filter if desired
        noiseIntervalResponse = filterAndVectorizeResponse(p.Results.filterMethod,noiseIntervalResponse,learningStructure);
        signalIntervalResponse = filterAndVectorizeResponse(p.Results.filterMethod,signalIntervalResponse,learningStructure);
        
        % Simulate trial
        response = simulateTrial(p.Results.learningMethod,false,noiseIntervalResponse,signalIntervalResponse,learningStructure);
        
        % Decide whether response was correct and keep track
        if (response)
            nCorrects(cc) = nCorrects(cc)+1;
        end
    end
end

%% Fit and plot psychometric function
fitContrasts = logspace(log10(min(testContrasts)),log10(max(testContrasts)),100)';
[threshold,fitFractionCorrect,paramsValues] = singleThresholdExtraction(testContrasts,nCorrects/p.Results.nTrialsPerTestContrast,thresholdParams.criterionFraction,p.Results.nTrialsPerTestContrast,fitContrasts);
hFig = figure; hold on
set(gca,'FontSize',rParams.plotParams.axisFontSize);
plot(log10(testContrasts), nCorrects/p.Results.nTrialsPerTestContrast,'ro', 'MarkerSize', rParams.plotParams.markerSize, 'MarkerFaceColor', [1.0 0.5 0.50]);
plot(log10(fitContrasts),fitFractionCorrect,'r','LineWidth', 2.0);
plot(log10(threshold)*[1 1],[0 thresholdParams.criterionFraction],'b', 'LineWidth', 2.0);
axis 'square'
set(gca, 'YLim', [0 1.0],'XLim', log10([testContrasts(1) testContrasts(end)]), 'FontSize', 14);
xlabel('contrast', 'FontSize' ,rParams.plotParams.labelFontSize, 'FontWeight', 'bold');
ylabel('percent correct', 'FontSize' ,rParams.plotParams.labelFontSize, 'FontWeight', 'bold');
box off; grid on
title({sprintf('LMangle = %2.1f deg, LMthreshold (%0.4f%%,%0.4f%%)', atan2(testConeContrasts(2), testConeContrasts(1))/pi*180, ...
    100*threshold*testConeContrasts(1), 100*threshold*testConeContrasts(2)) ; ''}, ...
    'FontSize',rParams.plotParams.titleFontSize);
%rwObject.write(sprintf('LMPsychoFunctions_%d',ii),hFig,paramsList,writeProgram,'Type','figure');

end

function [response,learningStructure] = simulateTrial(method,learn,noiseIntervalResponse,signalIntervalResponse,learningStructure)
% [response,learningStructure] = simulateTrial(method,learn,noiseIntervalResponse,signalIntervalResponse,learningStructure)
%
% Simulate out a trial and update the learning structure, which is method
% dependent.

% Determine which trial type (noise-signal, signal-noise) this is,
% and construct trial response
trialType = round(rand);
if (trialType == 0)
    trialResponse = [noiseIntervalResponse ; signalIntervalResponse];
else
    trialResponse = [signalIntervalResponse ; noiseIntervalResponse];
end

% Classify trial response.
switch (method)
    
    case {'nnmse' 'nncorr'}
        % Nearest neighbor methodsJust guess if we haven't seen exemplars
        % of each type.  Use nearest neighbor classification on trials we
        % have seen if we have.
        if (learningStructure.nTrials0 == 0 & learningStructure.nTrials1 == 0)
            responseType = round(rand);
        else
            trialDistances0 = zeros(learningStructure.nTrials0,1);
            for ii = 1:learningStructure.nTrials0
                switch (method)
                    case 'nnmse'
                        trialDistances0(ii) = norm(trialResponse-learningStructure.trialResponseData0(:,ii));
                    case 'nncorr'
                        trialDistances0(ii)  = 1/corr(trialResponse,learningStructure.trialResponseData0(:,ii));
                end
            end
            trialDistances1 = zeros(learningStructure.nTrials1,1);
            for ii = 1:learningStructure.nTrials1
                switch (method)
                    case 'nnmse'
                        trialDistances1(ii) = norm(trialResponse-learningStructure.trialResponseData1(:,ii));
                    case 'nncorr'
                        trialDistances1(ii) = 1/corr(trialResponse,learningStructure.trialResponseData1(:,ii));
                end
            end
            if (min(trialDistances0) < min(trialDistances1))
                responseType = 0;
            else
                responseType = 1;
            end
        end
        
    case 'uniflikely'
        % Figure out which interval is more likely to be a uniform field
        % and respond accordingly.
        responseDim = length(trialResponse);
        firstIntervalResponse = trialResponse(1:responseDim/2);
        secondIntervalResponse = trialResponse(responseDim/2+1:responseDim);
        firstLikely = computeUnifLikelihood(firstIntervalResponse,learningStructure);
        secondLikely = computeUnifLikelihood(secondIntervalResponse,learningStructure);
        if (trialType == 1)
            if (secondLikely < firstLikely)
                responseType = 1;
            else
                responseType = 0;
            end
        else
            if (secondLikely < firstLikely)
                responseType = 0;
            else
                responseType = 1;
            end
        end
    otherwise
        error('Unknown learning method specified');
end

% Learn if desired
%
% Save trial response data for future classification
if (learn)
    if (trialType == 0)
        learningStructure.nTrials0 = learningStructure.nTrials0+1;
        learningStructure.trialResponseData0(:,learningStructure.nTrials0) = trialResponse;
    else
        learningStructure.nTrials1 = learningStructure.nTrials1+1;
        learningStructure.trialResponseData1(:,learningStructure.nTrials1) = trialResponse;
    end
end

% Decide whether response was correct or not
if (responseType == trialType)
    response = 1;
else
    response = 0;
end
end

function logLikely = computeUnifLikelihood(intervalResponse,learningStructure)
% logLikely = computeUnifLikelihood(intervalResponse,learningStructure)
%
% Compute the log likelihood that the respones from an interval came from a
% uniform field, based on a simple normal approximation.

logLikely = 0;
for ii = 1:3
    index = learningStructure.coneIndex{ii};
    u = mean(intervalResponse(index));
    s = std(intervalResponse(index));
    z = (intervalResponse(index)-u)/s;
    logLikely = logLikely + sum(log(normpdf(z)));
end
end


function filteredResponse = filterAndVectorizeResponse(method,response,learningStructure)
% function filteredResponse = filterResponse(method,response,learningStructure)
%
% Optionally apply a lowpass filter and downsampling to mosaic responses,
% and make it a column vector.

switch (method)
    case 'none'
        % Just vectorize from image format
        filteredResponse = response(:);
        
    case 'lowpass'
        % Lowpass filter and vectorize
        % A filter param is in learningStructure.lowpassParam (space constant in microns)
        % The mosaic object is in learningStructure.theMosaic
        % See t_coneMosaicLowPassResponses() for a detailed demonstration
        % of the lowPassMosaicResponse() coneMosaic method
        [filteredResponse, ~, ~, ~] = learningStructure.theMosaic.lowPassMosaicResponse(response, learningStructure.lowpassParam * [1 1 1]);
        
        % Vectorize from image format
        filteredResponse = filteredResponse(:);
end

end
