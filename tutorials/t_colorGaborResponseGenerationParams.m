function params = t_colorGaborRespnseGenerationParams(varargin)
% params = t_colorGaborRespnseGenerationParams(varargin)
%
% This tutorial specifies the parameter structure for generating
% responses in the form expected by other programs in this project.  It
% provides reasonable illustrative defaults, although in real calculations
% it would be natural to vary many of these parameters.

%% Type field for parameters
params.type = 'ResponseGeneration';
 
%% Define parameters of a spatial gabor pattern
%
%   fieldOfViewDegs - Field of view in degrees, horizontal direction.
%   cyclesPerDegree - Grating cycles per degree.
%   gaussianFWHMDegs - Full width at half max of spatial Gaussian window.
%   row - Row dimension of scene on monitor
%   col - Col dimension of scene on monitor
%   ang - Angle of grating, in radians
%   ph  - Phase of grating, in radians relative to image center
%   viewingDistance - Viewing distance of observer from monitor in meters.
params.gaborParams.type = 'Gabor';
params.gaborParams.fieldOfViewDegs = 4;
params.gaborParams.cyclesPerDegree = 2;
params.gaborParams.gaussianFWHMDegs = 1.5;
params.gaborParams.row = 128;
params.gaborParams.col = 128;
params.gaborParams.ang = 0;
params.gaborParams.ph = 0;
params.gaborParams.viewingDistance = 0.75;

%% Define color modulation parameters
%   contrast - Contrast specfied relative to coneContrasts.
%              Can be a vector of contrasts.
%   coneContrasts - Color direction of grating in cone contrast space
%                   Can be a 3 by N matrix of contrast directions.
%   backgroundxYY - Colorimetric specification of background, in CIE xyY (cd/m2)
%   monitorFile - Isetbio display description of monitor on which grating is shown.
%   leakageLum - Luminance when monitor input is zero.
params.colorModulationParams.type = 'ColorModulation';
params.colorModulationParams.contrast = 1;
params.colorModulationParams.coneContrasts = [0.05 -0.05 0]';
params.colorModulationParams.backgroundxyY = [0.27 0.30 49.8]';
params.colorModulationParams.monitorFile = 'CRT-MODEL';
gaborParams.leakageLum = 1.0;

%% Parameters related to temporal properties of stimulus and response
%
%   frameRate - Frame rate in Hz of display device.
%   windowTauInSeconds - Standard deviation of Gaussian stimulus window.
%   stimulusDurationInSeconds - Stimulus duration.
%   stimulusSamplingIntervalInSeconds - How often we sample the stimulus time sequence.
%   secondsToInclude - Portion of response movie to include for classification.
%   secondsToIncludeOffset - Temporal offset of included window.
%   eyeDoNotMove - Boolean, set to true for perfect fixation.
%   simulationTimeStepSecs - Time step used in temporal simulation
params.temporalParams.type = 'TemporalParams';
params.temporalParams.frameRate = 60;
params.temporalParams.windowTauInSeconds = 0.165;
params.temporalParams.stimulusDurationInSeconds = 5*params.temporalParams.windowTauInSeconds;
params.temporalParams.stimulusSamplingIntervalInSeconds = 1/params.temporalParams.frameRate;
params.temporalParams.secondsToInclude = 0.050;
params.temporalParams.secondsToIncludeOffset = 0;
params.temporalParams.eyesDoNotMove = false; 
params.temporalParams.simulationTimeStepSecs = 5/1000; 

% Optional CRT raster effects.
% 
% Some of our routines understand these.
%    addCRTrasterEffect  - Incorporate a simple model of CRT raster timing.
%    rasterSamples - Temporal samples per simulation time step.
%                    This is used to speed up the simulation.
params.temporalParams.addCRTrasterEffect = false;
params.temporalParams.rasterSamples = 5; 

% Some computed temporal parameters
[params.temporalParams.sampleTimes,params.temporalParams.gaussianTemporalWindow] = gaussianTemporalWindowCreate(params.temporalParams);
params.temporalParams.nSampleTimes = length(params.temporalParams.sampleTimes);
if (params.temporalParams.addCRTrasterEffect)
   params.simulationTimeStepSecs = params.simulationTimeStepSecs/params.temporalParams.rasterSamples;
end

%% Properties related to computing the retinal image
%
%  fieldOfViewDegs - Field of view computed
%  offAxis - Boolean, compute falloff of intensity with position
%  blur - Boolean, incorporate optical blurring
%  lens - Boolean, incorporate filtering by lens
params.oiParams.type = 'OpticsParams';
params.oiParams.fieldOfViewDegs = params.gaborParams.fieldOfViewDegs;
params.oiParams.offAxis = false;
params.oiParams.blur = true;
params.oiParams.lens = true;

%% Properties of the cone mosaic
%
%  fieldOfViewDegs - Field of view computed
%  macular - Boolean, include macular pigment (may not be implemeted yet)
%  LMSRatio - Vector giving L, M, S cone ratio
%  timeStepInSeconds - Time step to compute responses on
%  integrationTimeInSeconds - Cone integration time.  Generally the same as time step
%  isomerizationNoise - Boolean, add isomerization Poisson noise?
%  osModel - What model to use to compute photocurrent
params.mosaicParams.type = 'MosaicParams';
params.mosaicParams.fieldOfViewDegs = params.gaborParams.fieldOfViewDegs;
params.mosaicParams.macular = true;
params.mosaicParams.LMSRatio = [0.62 0.31 0.07];
params.mosaicParams.timeStepInSeconds = params.temporalParams.stimulusSamplingIntervalInSeconds;
params.mosaicParams.integrationTimeInSeconds = params.mosaicParams.timeStepInSeconds;
params.mosaicParams.isomerizationNoise = false;
params.mosaicParams.osModel = 'Linear';

%% Parameters for plots
params.plotParams.type = 'PlotParams';
params.plotParams.labelFontSize = 12;
params.plotParams.titleFontSize = 14;
params.plotParams.axisFontSize = 8;