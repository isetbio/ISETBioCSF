function dirname = paramsToColorModulationDirName(obj,params)
% pdirname = paramsToColorModulationDirName(obj,params)
% 
% Generate a directory names that captures the color stimulus
% root parameters, but skips the color direction and contrast info,
% because these have been summarized over.

if (~strcmp(params.type,'ColorModulationRoot'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('Lum%0.2f_lf%0.1f',...
    params.backgroundxyY(3),...
    params.lumFactor);
