function [sampleTimes, squareTemporalWindow, rasterModulation] = squareTemporalWindowCreate(temporalParams)
%
% Create a square temporal window.

if (temporalParams.stimulusDurationInSeconds == temporalParams.stimulusSamplingIntervalInSeconds)
    sampleTimes = [0];
    squareTemporalWindow = 1;
    rasterModulation = 1;
    return;
end


stimulusSamples = round(temporalParams.stimulusDurationInSeconds/temporalParams.stimulusSamplingIntervalInSeconds);
if (isfield(temporalParams, 'secondsForResponseStabilization'))
    stabilizingTimeSamples = round(temporalParams.secondsForResponseStabilization/temporalParams.stimulusSamplingIntervalInSeconds);
else
    stabilizingTimeSamples = 0;
end

if (isfield(temporalParams, 'secondsForResponseExtinction'))
    extinctionTimeSamples = round(temporalParams.secondsForResponseExtinction/temporalParams.stimulusSamplingIntervalInSeconds);
else
    extinctionTimeSamples = 0;
end

sampleTimes = (1:(stabilizingTimeSamples + stimulusSamples + extinctionTimeSamples))*temporalParams.stimulusSamplingIntervalInSeconds;
sampleTimes = sampleTimes - mean(sampleTimes);
sampleTimes = sampleTimes - (sampleTimes(2)-sampleTimes(1))/2;

squareTemporalWindow = zeros(1,numel(sampleTimes));
squareTemporalWindow(abs(sampleTimes) <= temporalParams.stimulusDurationInSeconds/2) = 1;
rasterModulation = [];
end

