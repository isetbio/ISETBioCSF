function temporalParams = temporalParamsGenerate(varargin)
% temporalParams = temporalParamsGenerate(varargin)
%
% Parameters related to temporal properties of stimulus and response
%
%   frameRate - Frame rate in Hz of display device.
%   windowTauInSeconds - Standard deviation of Gaussian stimulus window.
%   stimulusDurationInSeconds - Stimulus duration.
%   stimulusSamplingIntervalInSeconds - How often we sample the stimulus time sequence.
%   secondsToInclude - Portion of response movie to include for classification.
%   secondsToIncludeOffset - Temporal offset of included window.
%   eyeDoNotMove - Boolean, set to true for perfect fixation.
%   simulationTimeStepSecs - Time step used in temporal simulation

temporalParams.type = 'TemporalParams';

temporalParams.frameRate = 60;
temporalParams.windowTauInSeconds = 0.165;
temporalParams.stimulusDurationInSeconds = 5*temporalParams.windowTauInSeconds;
temporalParams.stimulusSamplingIntervalInSeconds = 1/temporalParams.frameRate;
temporalParams.secondsToInclude = 0.050;
temporalParams.secondsToIncludeOffset = 0;
temporalParams.eyesDoNotMove = false; 
temporalParams.simulationTimeStepSecs = 5/1000; 

% Optional CRT raster effects.
% 
% Some of our routines understand these.
%    addCRTrasterEffect  - Incorporate a simple model of CRT raster timing.
%    rasterSamples - Temporal samples per simulation time step.
%                    This is used to speed up the simulation.
temporalParams.addCRTrasterEffect = false;
temporalParams.rasterSamples = 5; 

% Some computed temporal parameters
[temporalParams.sampleTimes,temporalParams.gaussianTemporalWindow] = gaussianTemporalWindowCreate(temporalParams);
temporalParams.nSampleTimes = length(temporalParams.sampleTimes);

% if (temporalParams.addCRTrasterEffect)
%    simulationTimeStepSecs = simulationTimeStepSecs/temporalParams.rasterSamples;
% end