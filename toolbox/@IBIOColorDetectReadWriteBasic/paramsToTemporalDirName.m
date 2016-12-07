function dirname = paramsToTemporalDirName(obj,temporalParams)
% dirname = paramsToTemporalDirName(obj,temporalParams)
% 
% Generate a directory names that captures the temporal parameters used to
% generate the responses.

if (~strcmp(temporalParams.type,'Temporal')) && (~strcmp(temporalParams.type,'Temporal_v2'))
    error('Incorrect parameter type passed');
end

if (strcmp(temporalParams.type,'Temporal'))
    
    if (isfield(temporalParams, 'emPathType'))
        
    switch (temporalParams.emPathType)
        case 'Zero'
            nemNumber = 1;
        case 'Dynamic'
            nemNumber = 0;
        case 'Frozen'
            nemNumber = 2;
    end
    else
        nemNumber = temporalParams.eyesDoNotMove;
    end
    
    dirname = sprintf('tau%0.3f_dur%0.2f_nem%d_use%0.0f_off%0.0f',...
        temporalParams.windowTauInSeconds, ...
        temporalParams.stimulusDurationInSeconds, ...
        nemNumber, ...
        temporalParams.secondsToInclude, ...
        temporalParams.secondsToIncludeOffset);
else
    dirname = sprintf('[STIM_TEMPORAL]_rampDur%0.3f_rampTau%0.3f',...
        temporalParams.rampDurationSecs, ...
        temporalParams.rampTauSecs ...
        );
end

