function dirname = paramsToColorModulationDirName(obj,params)
% dirname = paramsToColorModulationDirName(obj,params)
%
% Generate a directory names that captures the color stimulus
% parameters

if (~strcmp(params.type,'ColorModulation')) && (~strcmp(params.type,'ColorModulation_v2'))
    error('Incorrect parameter type passed');
end

if (strcmp(params.type,'ColorModulation'))
    switch (params.modulationType)
        case 'monitor'
            dirname = sprintf('L%0.0f_M%0.0f_S%0.0f_con%0.5f',...
                100*params.coneContrasts(1),...
                100*params.coneContrasts(2),...
                100*params.coneContrasts(3),...
                100*params.contrast);
        case 'AO'
            dirname = sprintf('%d_%0.1f_con%0.5f',params.spotWavelengthNm,params.spotCornealPowerUW,100*params.contrast);
        otherwise
            error('Unknown background type specified');
    end
else
     switch (params.device)
        case 'Monitor'
            dirname = sprintf('[CONE_MODULATION]_unitVector%0.3f_%0.3f_%0.3f_strength%0.5f%',...
                params.coneContrastUnitVector(1),...
                params.coneContrastUnitVector(2),...
                params.coneContrastUnitVector(3),...
                params.stimulusStrength);
     otherwise
            error('Unknown background type specified');
    end
            
end




