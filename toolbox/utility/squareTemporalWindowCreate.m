function [sampleTimes, squareTemporalWindow, rasterModulation] = squareTemporalWindowCreate(temporalParams)
%
% Create a square temporal window.

if (temporalParams.stimulusDurationInSeconds == temporalParams.stimulusSamplingIntervalInSeconds)
    sampleTimes = [0];
    squareTemporalWindow = 1;
    rasterModulation = 1;
    return;
end

nPositiveTimeSamples = ceil(0.5*temporalParams.stimulusDurationInSeconds/temporalParams.stimulusSamplingIntervalInSeconds);
if (isfield(temporalParams, 'secondsForResponseStabilization'))
    stabilizingTimeSamples = ceil(0.5*temporalParams.secondsForResponseStabilization/temporalParams.stimulusSamplingIntervalInSeconds);
else
    stabilizingTimeSamples = 0;
end

if (isfield(temporalParams, 'secondsForResponseExtinction'))
    extinctionTimeSamples = ceil(0.5*temporalParams.secondsForResponseExtinction/temporalParams.stimulusSamplingIntervalInSeconds);
else
    extinctionTimeSamples = 0;
end

sampleTimes = linspace(...
    -(nPositiveTimeSamples+stabilizingTimeSamples)*temporalParams.stimulusSamplingIntervalInSeconds, ...
    (nPositiveTimeSamples+extinctionTimeSamples)*temporalParams.stimulusSamplingIntervalInSeconds, ...
    2*nPositiveTimeSamples+1+stabilizingTimeSamples+extinctionTimeSamples);

squareTemporalWindow = zeros(1,numel(sampleTimes));
squareTemporalWindow(abs(sampleTimes) <= temporalParams.stimulusDurationInSeconds/2) = 1;

sampleTimes = linspace(sampleTimes(1), sampleTimes(end), numel(squareTemporalWindow));
rasterModulation = [];

end

