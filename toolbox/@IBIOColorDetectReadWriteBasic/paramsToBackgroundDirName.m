function dirname = paramsToBackgroundDirName(obj,params)
% dirname = paramsToColorModulationDirName(obj,params)
%
% Generate a directory names that captures the background
% parameters

if (~strcmp(params.type,'Background')) && (~strcmp(params.type,'Background_v2'))
    error('Incorrect parameter type passed');
end

if (strcmp(params.type,'Background')) 
    switch (params.backgroundType)
        case 'monitor'
            dirname = sprintf('lum%0.2f_lf%0.1f_leak%0.1f',...
                params.backgroundxyY(3),...
                params.lumFactor,...
                params.leakageLum);
        case 'AO'
            for ii = 1:length(params.backgroundWavelengthsNm)
                dirnameTmp = sprintf('%d_%0.1f',params.backgroundWavelengthsNm(ii),params.backgroundCornealPowerUW(ii));
                if (ii == 1)
                    dirname = dirnameTmp;
                else
                    dirname = [dirname '_' dirnameTmp];
                end
            end
        otherwise
            error('Unknown background type specified');
    end
else
    switch (params.device)
        case 'Monitor'
            dirname = sprintf('[BACKGROUND]_Lum%0.2f',...
                params.backgroundxyY(3) ...
            );
    otherwise
            error('Unknown background type specified');        
    end
end






