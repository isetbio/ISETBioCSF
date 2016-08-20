function c_EffectOfTrainingSize(varargin)
% c_EffectOfTrainingSize(varargin)
%
% Study effect of SVM training size on estimated thresholds, and compare to
% signal known exactly ML classifier and empirical ML classifier.
%
% This looks at L+M detection thrsholds at a moderate spatial frequency (10
% cpd).
%
% Key/value pairs
%   'PlotOnly' - true/false (default false).  Just plot the summary figure.

%% Parse input
p = inputParser;
p.addParameter('PlotOnly',false,@islogical);
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
rParams.gaborParams.cyclesPerDegree = 10;
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
testDirectionParams.contrastScale = 'log';    % choose between 'linear' and 'log'

%% Parameters related to how we find thresholds from responses
%
% Use default
thresholdParams = thresholdParamsGenerate;

%% Loop over triaing samples
effectOfTrainingSize.nTrainingSamplesList = [50 100 500 1000 5000 10000];
if (~p.Results.PlotOnly)
    for tt = 1:length(effectOfTrainingSize.nTrainingSamplesList)
        
        % Set number of trials
        testDirectionParams.trialsNum = effectOfTrainingSize.nTrainingSamplesList(tt);
        
        %% Compute response instances
        t_colorGaborConeCurrentEyeMovementsResponseInstances('rParams',rParams,'testDirectionParams',testDirectionParams,'compute',true,'visualizeResponses',false);
        
        %% Find thresholds and summarize, template max likeli
        thresholdParams.method = 'mlpt';
        t_colorGaborDetectFindPerformance('rParams',rParams,'testDirectionParams',testDirectionParams,'thresholdParams',thresholdParams,'compute',true,'plotSvmBoundary',false,'plotPsychometric',false);
        effectOfTrainingSize.mlptThresholds(tt) = t_plotGaborDetectThresholdsOnLMPlane('rParams',rParams,'LMPlaneInstanceParams',testDirectionParams,'thresholdParams',thresholdParams,'plotPsychometric',false,'plotEllipse',false);
        close all;
        
        %% Find thresholds and summarize, svm
        thresholdParams.method = 'svm';
        thresholdParams.PCAComponents = 500;
        t_colorGaborDetectFindPerformance('rParams',rParams,'testDirectionParams',testDirectionParams,'thresholdParams',thresholdParams,'compute',true,'plotSvmBoundary',false,'plotPsychometric',false);
        effectOfTrainingSize.svmThresholds(tt) = t_plotGaborDetectThresholdsOnLMPlane('rParams',rParams,'LMPlaneInstanceParams',testDirectionParams,'thresholdParams',thresholdParams,'plotPsychometric',false,'plotEllipse',false);
        close all;
        
        %% Find thresholds and summarize, empirical max likeli
        thresholdParams.method = 'mlpe';
        t_colorGaborDetectFindPerformance('rParams',rParams,'testDirectionParams',testDirectionParams,'thresholdParams',thresholdParams,'compute',true,'plotSvmBoundary',false,'plotPsychometric',false);
        effectOfTrainingSize.mlpeThresholds(tt) = t_plotGaborDetectThresholdsOnLMPlane('rParams',rParams,'LMPlaneInstanceParams',testDirectionParams,'thresholdParams',thresholdParams,'plotPsychometric',false,'plotEllipse',false);
        close all;
    end
    
    %% Write out the data
    %
    % Set trialsNum to 0 to define a summary directory name
    fprintf('Writing performance data ... ');
    testDirectionParams.trialsNum = 0;
    paramsList = {rParams.gaborParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    rwObject.write('effectOfTrainingSize',effectOfTrainingSize,paramsList,writeProgram);
    fprintf('done\n');
else
    % Read in the data
    fprintf('Reading performance data ...');
    testDirectionParams.trialsNum = 0;
    paramsList = {rParams.gaborParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    effectOfTrainingSize = rwObject.read('effectOfTrainingSize',paramsList,writeProgram);
    fprintf('done\n');
end

%% Make a plot of estimated threshold versus training set size
%
% The way the plot is coded counts on the test contrasts never changing
% across the conditions, which we could explicitly check for here.
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
