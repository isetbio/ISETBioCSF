function dirname = paramsToLMSsamplingDirName(obj,params)
% dirname = paramsToLMSsamplingDirName(obj,params)
% 
% Generate a directory names that captures the session parameters.


if (~strcmp(params.type,'LMSsampling'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('[LMS_SPACE_SAMPLING]_azimuth%2.0f_elevation%2.0f_stimStrengthAxis_%0.5f_%0.5f_%0.0f',...
                params.azimuthAngle, ...
                params.elevationAngle, ...
                params.stimulusStrengthAxis(1), ...
                params.stimulusStrengthAxis(end), ...
                numel(params.stimulusStrengthAxis) ...
                );
end




