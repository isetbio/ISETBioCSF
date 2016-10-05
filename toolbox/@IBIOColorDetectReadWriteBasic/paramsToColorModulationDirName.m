function dirname = paramsToColorModulationDirName(obj,params)
% dirname = paramsToColorModulationDirName(obj,params)
%
% Generate a directory names that captures the color stimulus
% parameters

if (~strcmp(params.type,'ColorModulation'))
    error('Incorrect parameter type passed');
end

switch (params.modulationType)
    case 'monitor'
        dirname = sprintf('L%0.0f_M%0.0f_S%0.0f_con%0.5f',...
            100*params.coneContrasts(1),...
            100*params.coneContrasts(2),...
            100*params.coneContrasts(3),...
            100*params.contrast);
    case 'AO'
        dirname = sprintf('%d_%0.1f_con%0.5f',params.spotWavelengthNm,params.spotCornealIrradianceUW,100*params.contrast);
    otherwise
        error('Unknown background type specified');
end



