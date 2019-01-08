function [sampleTimes, squareTemporalWindow, rasterModulation] = squareTemporalWindowCreate(temporalParams)
%
% Special case where we generate a temporal window with a singe time sample
if (temporalParams.stimulusDurationInSeconds == temporalParams.stimulusSamplingIntervalInSeconds) && ...
   (temporalParams.singleBinTemporalWindowIfPossible)
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

sampleTimes = 1:(stabilizingTimeSamples + stimulusSamples + extinctionTimeSamples);
sampleTimes = sampleTimes - stabilizingTimeSamples - round(stimulusSamples/2);
sampleTimes = sampleTimes * temporalParams.stimulusSamplingIntervalInSeconds;

squareTemporalWindow = zeros(1,numel(sampleTimes));
onTime = find(sampleTimes >= -temporalParams.stimulusDurationInSeconds/2);
squareTemporalWindow(onTime(1)+(0:(stimulusSamples-1))) = 1;

rasterModulation = [];
end


% BELOW BETTER VERSIOn BUT NEEDS TESTING
function [sampleTimes, squareTemporalWindow, rasterModulation] = squareTemporalWindowCreateUpdate(temporalParams)
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

totalSamples = stabilizingTimeSamples + stimulusSamples + extinctionTimeSamples
totalTime = totalSamples*temporalParams.stimulusSamplingIntervalInSeconds

sampleTimes = 1:totalSamples;
sampleTimes = sampleTimes - round(totalSamples/2)-0.5;
sampleTimes = sampleTimes * temporalParams.stimulusSamplingIntervalInSeconds;

squareTemporalWindow = zeros(1,numel(sampleTimes));
onTime = find(sampleTimes > - (temporalParams.stimulusDurationInSeconds/2+temporalParams.stimulusSamplingIntervalInSeconds/2));
squareTemporalWindow(onTime(1)+(0:(stimulusSamples-1))) = 1;

rasterModulation = [];
end

