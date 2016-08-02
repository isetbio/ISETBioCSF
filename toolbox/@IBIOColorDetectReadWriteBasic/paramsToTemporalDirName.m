function dirname = paramsToResponseGenerationDirName(obj,temporalParams)
% dirname = paramsToResponseGenerationDirName(obj,temporalParams)
% 
% Generate a directory names that captures the temporal parameters used to
% generate the responses.

if (~strcmp(temporalParams.type,'Temporal'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('tau%0.3f_dur%0.2f_nem%0.0f_use%0.0f_off%0.0f',...
    temporalParams.windowTauInSeconds, ...
    temporalParams.stimulusDurationInSeconds, ...
    temporalParams.eyesDoNotMove, ...
    temporalParams.secondsToInclude, ...
    temporalParams.secondsToIncludeOffset);

