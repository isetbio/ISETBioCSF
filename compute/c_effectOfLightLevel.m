function c_effectOfLightLevel
% c_effectOfLightLevel
%
% Compute threshold ellipses at multiple overall light levels.  We should
% see contrast thresholds decrease at the isomerizations, but remain more
% or less constant at the photocurrent.

%% Clear
ieInit; close all;

%% Get the parameters we need
%
% Start with default
rParams = colorGaborResponseParamsGenerate;

% Override some defult parameters
rParams.temporalParams.simulationTimeStepSecs = 10/1000;
rParams.temporalParams.stimulusDurationInSeconds = 1;
rParams.temporalParams.stimulusSamplingIntervalInSeconds = 1/60;
rParams.temporalParams.secondsToInclude = rParams.temporalParams.stimulusDurationInSeconds;
rParams.temporalParams.eyesDoNotMove = true;

rParams.mosaicParams.timeStepInSeconds = rParams.temporalParams.simulationTimeStepSecs;
rParams.mosaicParams.integrationTimeInSeconds = rParams.mosaicParams.timeStepInSeconds;
rParams.mosaicParams.isomerizationNoise = true;
rParams.mosaicParams.osNoise = true;
rParams.mosaicParams.osModel = 'Linear';

%% Parameters that define the LM instances we'll generate here
%
% Use default LMPlane.
testDirectionParams = LMPlaneInstanceParamsGenerate;

%% Parameters related to how we find thresholds from responses
%
% Use default
thresholdParams = thresholdParamsGenerate;

%% Compute response instances
t_colorGaborConeCurrentEyeMovementsResponseInstances(rParams,testDirectionParams);

%% Find thresholds
t_colorGaborDetectFindPerformance(rParams,testDirectionParams,thresholdParams);

%% LMPlane summary
t_plotGaborDetectThresholdsOnLMPlane(rParams,LMPlaneInstanceParams,thresholdParams);
