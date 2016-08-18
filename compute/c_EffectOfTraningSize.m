function c_EffectOfTrainingSize
% c_EffectOfTrainingSize
%
% Study effect of SVM training size on estimated thresholds, and compare to
% signal known exactly ML classifier and empirical ML classifier.
%
% This looks at L+M detection thrsholds at a moderate spatial frequency (10
% cpd).

%% Clear
ieInit; close all;

%% Loop over triaing samples
%effectOfTrainingSize.nTrainingSamplesList = [50 100 500 1000 5000];
effectOfTrainingSize.nTrainingSamplesList = [5000];
for tt = 1:length(effectOfTrainingSize.nTrainingSamplesList)
    
    %% Get the parameters we need
    %
    % Start with default
    rParams = colorGaborResponseParamsGenerate;
    
    % Override some defult parameters
    
    % Get stimulus parameters correct
    %
    % The stimulus was half-cosine windowed to contain 7.5 cycles.  We set
    % our half-cosine window to match that and also make the field of view
    % just a tad bigger.
    rParams.gaborParams.windowType = 'halfcos';
    rParams.gaborParams.cyclesPerDegree = 10;
    rParams.gaborParams.gaussianFWHMDegs = 3.75*(1/rParams.gaborParams.cyclesPerDegree);
    rParams.gaborParams.fieldOfViewDegs = 2.1*rParams.gaborParams.gaussianFWHMDegs;
    
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
    testDirectionParams.trialsNum = effectOfTrainingSize.nTrainingSamplesList(tt);
    
    % Number of contrasts to run in each color direction
    testDirectionParams.nContrastsPerDirection = 20;
    testDirectionParams.lowContrast = 0.0001;
    testDirectionParams.highContrast = 0.1;
    testDirectionParams.contrastScale = 'log';    % choose between 'linear' and 'log'
    
    %% Parameters related to how we find thresholds from responses
    %
    % Use default
    thresholdParams = thresholdParamsGenerate;
    
    %% Compute response instances
%     t_colorGaborConeCurrentEyeMovementsResponseInstances(rParams,testDirectionParams);
    
    %% Find thresholds and summarize, template max likeli
%     thresholdParams.method = 'mlpt';
%     t_colorGaborDetectFindPerformance(rParams,testDirectionParams,thresholdParams);
%     effectOfTrainingSize.mlptThresholds(tt) = t_plotGaborDetectThresholdsOnLMPlane(rParams,testDirectionParams,thresholdParams);
%     close all;
    
    %% Find thresholds and summarize, svm
    thresholdParams.method = 'svm';
    thresholdParams.PCAComponents = 500;
    t_colorGaborDetectFindPerformance(rParams,testDirectionParams,thresholdParams);
    effectOfTrainingSize.svmThresholds(tt) = t_plotGaborDetectThresholdsOnLMPlane(rParams,testDirectionParams,thresholdParams);
    close all;
    
    %% Find thresholds and summarize, empirical max likeli
    thresholdParams.method = 'mlpe';
    t_colorGaborDetectFindPerformance(rParams,testDirectionParams,thresholdParams);
    effectOfTrainingSize.mlpeThresholds(tt) = t_plotGaborDetectThresholdsOnLMPlane(rParams,testDirectionParams,thresholdParams);
    close all;
end

% %% Write out the data
% fprintf('Writing performance data ... ');
% rwObject = IBIOColorDetectReadWriteBasic;
% writeProgram = mfilename;
% paramsList = {rParams.gaborParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
% rwObject.write('effectOfTrainingSize',effectOfTrainingSize,paramsList,writeProgram);
% fprintf('done\n');
% 
% %% Make a plot of estimated threshold versus training set size
% %
% % The way the plot is coded counts on the test contrasts never changing
% % across the conditions, which we could explicitly check for here.
% hFig = figure; clf; hold on
% set(gca,'FontSize', rParams.plotParams.axisFontSize);
% plot(effectOfTrainingSize.nTrainingSamplesList,[effectOfTrainingSize.mlptThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1),'r');
% plot(effectOfTrainingSize.nTrainingSamplesList,[effectOfTrainingSize.mlpeThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1),'g');
% plot(effectOfTrainingSize.nTrainingSamplesList,[effectOfTrainingSize.svmThresholds(:).thresholdContrasts]*effectOfTrainingSize.mlptThresholds(1).testConeContrasts(1),'b');
% xlabel('Training Set Size', 'FontSize' ,rParams.plotParams.labelFontSize, 'FontWeight', 'bold');
% ylabel('Threshold Contrast', 'FontSize' ,rParams.plotParams.labelFontSize, 'FontWeight', 'bold');
% box off; grid on
% title('Effect of Training Set Size','FontSize',rParams.plotParams.titleFontSize);
% rwObject.write('effectOfTrainingSize',hFig,paramsList,writeProgram,'Type','figure');
