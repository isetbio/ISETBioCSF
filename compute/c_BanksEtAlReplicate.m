function c_effectOfLightLevel
% c_effectOfLightLevel
%
% Compute threshold ellipses to replicate Banks et al, 1987,
% more or less.

%% Clear
ieInit; close all;

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
rParams.gaborParams.cyclesPerDegree = 3;
rParams.gaborParams.gaussianFWHMDegs = 3.25*(1/rParams.gaborParams.cyclesPerDegree);
rParams.gaborParams.fieldOfViewDegs = 2.1*rParams.gaborParams.gaussianFWHMDegs;

% Set background luminance
baseLum = 50;
backgroundParams.backgroundxyY = [0.33 0.33 baseLum]';
backgroundParams.monitorFile = 'CRT-MODEL';
backgroundParams.leakageLum = 1.0;
backgroundParams.lumFactor = 304/baseLum;

% Pupil size.  They used a 2mm artificial pupil
oiParams.pupilDiamMm = 2;

% Set duration equal to sampling interval to do just one frame.
%
% Their intervals were 100 msec each.
rParams.temporalParams.simulationTimeStepSecs = 100/1000;
rParams.temporalParams.stimulusDurationInSeconds = rParams.temporalParams.simulationTimeStepSecs;
rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.simulationTimeStepSecs;
rParams.temporalParams.secondsToInclude = rParams.temporalParams.simulationTimeStepSecs;
rParams.temporalParams.eyesDoNotMove = true;

% Their main calculation was without eye movements
temporalParams.eyesDoNotMove = true; 

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
testDirectionParams.deltaAngle = 45;

% Number of contrasts to run in each color direction
params.nContrastsPerDirection = 10; 
params.lowContrast = 0.001;
params.highContrast = 0.15;
params.contrastScale = 'log';    % choose between 'linear' and 'log'

%% Parameters related to how we find thresholds from responses
%
% Use default
thresholdParams = thresholdParamsGenerate;

%% Compute response instances
%t_colorGaborConeCurrentEyeMovementsResponseInstances(rParams,testDirectionParams);

%% Find thresholds
t_colorGaborDetectFindPerformance(rParams,testDirectionParams,thresholdParams);

%% LMPlane summary
t_plotGaborDetectThresholdsOnLMPlane(rParams,LMPlaneInstanceParams,thresholdParams);
