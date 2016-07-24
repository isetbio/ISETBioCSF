function dirname = paramsToColorModulationDirName(obj,params)
% pdirname = paramsToColorModulationDirName(obj,params)
% 
% Generate a directory names that captures the color stimulus
% parameters

if (~strcmp(params.type,'ColorModulation'))
    error('Incorrect parameter type passed');
end

theColorModulationName = sprintf('L%0.0f_M%0.0f_S%0.0f_con%0.1f_lum%0.2f',...
    100*params.coneContrasts(1),...
    100*params.coneContrasts(2),...
    100*params.coneContrasts(3),...
    100*params.contrast,...
    params.backgroundxyY(3));

dirname = [theColorModulationName];
