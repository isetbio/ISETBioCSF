function dirname = paramsToResponseSubsamplingDirName(obj,params)
% dirname = paramsToResponseSubsamplingDirName(obj,params)
% 
% Generate a directory names that captures the response subsampling params

if (~strcmp(params.type,'ResponseSubsampling'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('[RESPONSE_SUBSAMPLING]_durationMillisecs%0.1f_offsetMillisecs_%0.1f', ...
                1000.0*params.secondsToInclude, ...
                1000.0*params.secondsToIncludeOffset ...
          );

end




