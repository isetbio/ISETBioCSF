function c_EffectOfTrainingSize(varargin)
% c_EffectOfTrainingSize(varargin)
%
% Study effect of SVM training size on estimated thresholds, and compare to
% signal known exactly ML classifier and empirical ML classifier.
%
% This looks at L+M detection thrsholds, by default at a moderate spatial frequency (10
% cpd).  The stimulus size is inversely proportional to spatial frequency,
% so using a larger spatial frequency speeds things up.
%
% Key/value pairs
%   'nTrainingSamplesList' - vector (default [50 100 500 1000 5000].  Vector
%       of number of training samples to cycle through.
%   'cyclesPerDegree' - value (default 10). Spatial frequency of grating to be investigated.
%       This is reciprocally related to the size of the image/mosaic.
%   'computeResponses' - true/false (default true).  Compute responses.
%   'computeMLPTPerformance' - true/false (default true).  Compute template
%       maximum likelihood performance.
%   'computeMLPEPerformance' - true/false (default true).  Compute
%       empirical maximum likelihood performance.
%   'computeSVMPerformance' - true/false (default true).  Compute SVM performance.
%   'nPCAComponents' - integer (default 500).  Number of PCA components to
%       use with SVM.
%   'computeMLPTThresholds' - true/false (default true).  Compute template
%       maximum likelihood thresholds.
%   'computeMLPEThresholds' - true/false (default true).  Compute
%       empirical maximum likelihood thresholds.cl
%   'computeSVMThresholds' - true/false (default true).  Compute SVM thresholds.
%   'plotPsychometric' - true/false (default true).  Plot psychometric
%       functions.
%   'plotTrainingSize' - true/false (default true).  Plot results.

%% Parse input
p = inputParser;
p.addParameter('nTrainingSamplesList',[50 100 500 1000 5000],@isnumeric);
p.addParameter('cyclesPerDegree',10,@isnumeric);
p.addParameter('computeResponses',true,@islogical);
p.addParameter('computeMLPTPerformance',true,@islogical);
p.addParameter('computeMLPEPerformance',true,@islogical);
p.addParameter('computeSVMPerformance',true,@islogical);
p.addParameter('nPCAComponents',500,@isnumeric);
p.addParameter('computeMLPTThresholds',true,@islogical);
p.addParameter('computeMLPEThresholds',true,@islogical);
p.addParameter('computeSVMThresholds',true,@islogical);
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
testDirectionParams.contrastScale = 'log';    % choose between 'linear' and 'log'

%% Parameters related to how we find thresholds from responses
%
% Use default
thresholdParams = thresholdParamsGenerate;

%% Loop over triaing samples
effectOfTrainingSize.nTrainingSamplesList = p.Results.nTrainingSamplesList;
for tt = 1:length(effectOfTrainingSize.nTrainingSamplesList)
    
    % Set number of trials
    testDirectionParams.trialsNum = effectOfTrainingSize.nTrainingSamplesList(tt);
    
    %% Compute response instances
    if (p.Results.computeResponses)
        t_colorGaborConeCurrentEyeMovementsResponseInstances('rParams',rParams,'testDirectionParams',testDirectionParams,'compute',true,'visualizeResponses',false);
    end
    
    %% Find performance, template max likeli
    if (p.Results.computeMLPTPerformance)
        thresholdParams.method = 'mlpt';
        t_colorGaborDetectFindPerformance('rParams',rParams,'testDirectionParams',testDirectionParams,'thresholdParams',thresholdParams,'compute',true,'plotSvmBoundary',false,'plotPsychometric',false);
    end
    
     %% Find performance\, empirical max likeli
    if (p.Results.computeMLPEPerformance)
        thresholdParams.method = 'mlpe';
        t_colorGaborDetectFindPerformance('rParams',rParams,'testDirectionParams',testDirectionParams,'thresholdParams',thresholdParams,'compute',true,'plotSvmBoundary',false,'plotPsychometric',false);
    end
    
    %% Find performance, svm
    if (p.Results.computeSVMPerformance)
        thresholdParams.method = 'svm';
        thresholdParams.PCAComponents = p.Results.nPCAComponents;
        t_colorGaborDetectFindPerformance('rParams',rParams,'testDirectionParams',testDirectionParams,'thresholdParams',thresholdParams,'compute',true,'plotSvmBoundary',false,'plotPsychometric',false);
    end
   
    %% Fit psychometric functions
    if (p.Results.fitPsychometric)
        if (p.Results.computeMLPTThresholds)
            thresholdParams.method = 'mlpt';
            effectOfTrainingSize.mlptThresholds(tt) = t_plotGaborDetectThresholdsOnLMPlane('rParams',rParams,'LMPlaneInstanceParams',testDirectionParams,'thresholdParams',thresholdParams, ...
                'plotPsychometric',p.Results.plotPsychometric,'plotEllipse',false);
            close all;
        end
        
        if (p.Results.computeMLPEThresholds)      
            thresholdParams.method = 'mlpe';
            effectOfTrainingSize.mlpeThresholds(tt) = t_plotGaborDetectThresholdsOnLMPlane('rParams',rParams,'LMPlaneInstanceParams',testDirectionParams,'thresholdParams',thresholdParams, ...
                'plotPsychometric',p.Results.plotPsychometric,'plotEllipse',false);
            close all;
        end
        
        if (p.Results.computeSVMThresholds) 
            thresholdParams.method = 'svm';
            thresholdParams.PCAComponents = p.Results.nPCAComponents;
            effectOfTrainingSize.svmThresholds(tt) = t_plotGaborDetectThresholdsOnLMPlane('rParams',rParams,'LMPlaneInstanceParams',testDirectionParams,'thresholdParams',thresholdParams, ...
                'plotPsychometric',p.Results.plotPsychometric,'plotEllipse',false);
            close all;
        end
    end
end

%% Write out the data
%
% Set trialsNum to 0 to define a summary directory name
fprintf('Writing performance data ... ');
testDirectionParams.trialsNum = 0;
paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
rwObject = IBIOColorDetectReadWriteBasic;
writeProgram = mfilename;
rwObject.write('effectOfTrainingSize',effectOfTrainingSize,paramsList,writeProgram);
fprintf('done\n');

%% Make a plot of estimated threshold versus training set size
%
% The way the plot is coded counts on the test contrasts never changing
% across the conditions, which we could explicitly check for here.
if (p.Results.plotTrainingSize)
    fprintf('Reading performance data ...');
    testDirectionParams.trialsNum = 0;
    paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    effectOfTrainingSize = rwObject.read('effectOfTrainingSize',paramsList,writeProgram);
    fprintf('done\n');
    
    hFig = figure; clf; hold on
    fontBump = 4;
    set(gca,'FontSize', rParams.plotParams.axisFontSize+fontBump);
    plot(log10(effectOfTrainingSize.nTrainingSamplesList),[effectOfTrainingSize.mlptThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1), ...
        'ro','MarkerSize',rParams.plotParams.markerSize,'MarkerFaceColor','r');
    plot(log10(effectOfTrainingSize.nTrainingSamplesList),[effectOfTrainingSize.mlpeThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1), ...
        'go','MarkerSize',rParams.plotParams.markerSize,'MarkerFaceColor','g');
    plot(log10(effectOfTrainingSize.nTrainingSamplesList),[effectOfTrainingSize.svmThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1), ...
        'bo','MarkerSize',rParams.plotParams.markerSize,'MarkerFaceColor','b');
    plot(log10(effectOfTrainingSize.nTrainingSamplesList),[effectOfTrainingSize.mlptThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1), ...
        'r','LineWidth',rParams.plotParams.lineWidth);
    plot(log10(effectOfTrainingSize.nTrainingSamplesList),[effectOfTrainingSize.mlpeThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1), ...
        'g','LineWidth',rParams.plotParams.lineWidth);
    plot(log10(effectOfTrainingSize.nTrainingSamplesList),[effectOfTrainingSize.svmThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1), ...
        'b','LineWidth',rParams.plotParams.lineWidth);
    xlabel('Log10 Training Set Size', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
    ylabel('Threshold Contrast', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
    xlim([1 4]); ylim([0 1e-2]);
    legend({'MaxLikelihood, signal known', 'MaxLikelihood, signal estimated', sprintf('SVM, PCA %d',thresholdParams.PCAComponents)},'Location','NorthEast','FontSize',rParams.plotParams.labelFontSize+fontBump);
    box off; grid on
    title(sprintf('Effect of Training Set Size, %0.2f deg mosaic',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
    rwObject.write('effectOfTrainingSize',hFig,paramsList,writeProgram,'Type','figure');
end
