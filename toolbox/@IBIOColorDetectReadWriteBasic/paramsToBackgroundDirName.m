function dirname = paramsToBackgroundDirName(obj,params)
% dirname = paramsToColorModulationDirName(obj,params)
% 
% Generate a directory names that captures the background
% parameters

if (~strcmp(params.type,'Background'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('lum%0.2f_lf%0.1f_leak%0.1f',...
    params.backgroundxyY(3),...
    params.lumFactor,...
    params.leakageLum);

