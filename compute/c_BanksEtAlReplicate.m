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
%
% We start with a base luminance that we know is about mid-gray on the
% monitor we specify.  To change luminance, we specify a scale factor.
% This is eventually applied both to the background luminance and to the
% monitor channel spectra, so that we don't get unintersting out of gamut errors.
baseLum = 50;
theLum = 304;
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
testDirectionParams.trialsNum = 1000;

% Number of contrasts to run in each color direction
testDirectionParams.nContrastsPerDirection = 10; 
testDirectionParams.lowContrast = 0.001;
testDirectionParams.highContrast = 0.15;
testDirectionParams.contrastScale = 'log';    % choose between 'linear' and 'log'

%% Parameters related to how we find thresholds from responses
%
% Use default
thresholdParams = thresholdParamsGenerate;
thresholdParams.PCAComponents = 200;

%% Compute response instances
t_colorGaborConeCurrentEyeMovementsResponseInstances(rParams,testDirectionParams);

%% Find thresholds
t_colorGaborDetectFindPerformance(rParams,testDirectionParams,thresholdParams);

%% LMPlane summary
t_plotGaborDetectThresholdsOnLMPlane(rParams,testDirectionParams,thresholdParams);
