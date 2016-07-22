function params = t_colorGaborRespnseGenerationParams(varargin)
% params = t_colorGaborRespnseGenerationParams(varargin)
%
% This tutorial specifies the parameter structure for generating
% responses in the form expected by other programs in this project.  It
% provides reasonable illustrative defaults, although in real calculations
% it would be natural to vary many of these parameters.
 
%% Define parameters of a color gabor pattern
%
%   fieldOfViewDegs - Field of view in degrees, horizontal direction.
%   cyclesPerDegree - Grating cycles per degree.
%   gaussianFWHMDegs - Full width at half max of spatial Gaussian window.
%   row - Row dimension of scene on monitor
%   col - Col dimension of scene on monitor
%   contrast - Contrast specfied relative to coneContrasts.
%   ang - Angle of grating, in radians
%   ph  - Phase of grating, in radians relative to image center
%   coneContrasts - Color direction of grating in cone contrast space
%   backgroundxYY - Colorimetric specification of background, in CIE xyY (cd/m2)
%   monitorFile - Isetbio display description of monitor on which grating is shown.
%   viewingDistance - Viewing distance of observer from monitor in meters.
params.gaborParams.fieldOfViewDegs = 4;
params.gaborParams.cyclesPerDegree = 2;
params.gaborParams.gaussianFWHMDegs = 1.5;
params.gaborParams.row = 128;
params.gaborParams.col = 128;
params.gaborParams.contrast = 1;
params.gaborParams.ang = 0;
params.gaborParams.ph = 0;
params.gaborParams.coneContrasts = [0.05 -0.05 0]';
params.gaborParams.backgroundxyY = [0.27 0.30 49.8]';
params.gaborParams.monitorFile = 'CRT-MODEL';
params.gaborParams.viewingDistance = 0.75;

%% Parameters related to temporal properties of stimulus and response
%
%   frameRate - Frame rate in Hz of display device.
%   windowTauInSeconds - Standard deviation of Gaussian stimulus window.
%   stimulusDurationInSeconds - Stimulus duration.
%   stimulusSamplingIntervalInSeconds - How often we sample the stimulus time sequence.
%   secondsToInclude - Portion of response movie to include for classification.
%   secondsToIncludeOffset - Temporal offset of included window.
%   eyeDoNotMove - Boolean, set to true for perfect fixation.
params.temporalParams.frameRate = 60;
params.temporalParams.windowTauInSeconds = 0.165;
params.temporalParams.stimulusDurationInSeconds = 5*params.temporalParams.windowTauInSeconds;
params.temporalParams.stimulusSamplingIntervalInSeconds = 1/params.temporalParams.frameRate;
params.temporalParams.secondsToInclude = 0.050;
params.temporalParams.secondsToIncludeOffset = 0;
params.temporalParams.eyesDoNotMove = false; 

% Some computed temporal parameters
[params.temporalParams.sampleTimes,params.temporalParams.gaussianTemporalWindow] = gaussianTemporalWindowCreate(params.temporalParams);
params.temporalParams.nSampleTimes = length(params.temporalParams.sampleTimes);

%% Properties related to computing the retinal image
%
%  fieldOfViewDegs - Field of view computed
%  offAxis - Boolean, compute falloff of intensity with position
%  blur - Boolean, incorporate optical blurring
%  lens - Boolean, incorporate filtering by lens
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
params.mosaicParams.fieldOfViewDegs = params.gaborParams.fieldOfViewDegs;
params.mosaicParams.macular = true;
params.mosaicParams.LMSRatio = [0.6 0.3 0.1];
params.mosaicParams.timeStepInSeconds = params.temporalParams.stimulusSamplingIntervalInSeconds;
params.mosaicParams.integrationTimeInSeconds = params.mosaicParams.timeStepInSeconds;
params.mosaicParams.isomerizationNoise = false;
params.mosaicParams.osModel = 'Linear';

%% Parameters for plots
params.plotParams.labelFontSize = 12;
params.plotParams.titleFontSize = 14;
params.plotParams.axisFontSize = 8;