function dirname = paramsToSessionDirName(obj,params)
% dirname = paramsToSessionDirName(obj,params)
% 
% Generate a directory names that captures the session parameters.

if (~strcmp(params.type,'Session'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('%s',...
                params.sessionType...
                );
end




