function dirname = paramsToTemporalDirName(obj,temporalParams)
% dirname = paramsToTemporalDirName(obj,temporalParams)
% 
% Generate a directory names that captures the temporal parameters used to
% generate the responses.

if (~strcmp(temporalParams.type,'Temporal'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('ST_rampTauMilliSecs%0.0f_stimDurMilliSecs%0.0f_analyzedResponseDurationMilliSecs%0.0f_withOffsetMilliSecs%0.0f_eyeMovementPath%s',...
    temporalParams.windowTauInSeconds*1000, ...
    temporalParams.stimulusDurationInSeconds*1000, ...
    temporalParams.secondsToInclude*1000, ...
    temporalParams.secondsToIncludeOffset*1000, ...
    sprintf('%s%s', upper(temporalParams.emPathType(1)), temporalParams.emPathType(2:end)));
end

