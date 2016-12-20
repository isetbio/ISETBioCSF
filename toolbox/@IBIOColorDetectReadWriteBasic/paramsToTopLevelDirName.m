function dirname = paramsToTopLevelDirName(obj,topLevelDirParams)
% dirname = paramsToTopLevelDirName(obj,topLevelDirParams)
% 
% Generate a directory name that captures the topLevelDir parameters.

if (~strcmp(topLevelDirParams.type,'TopLevelDir'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('[%s]', topLevelDirParams.name);

