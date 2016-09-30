function dirname = paramsToBackgroundDirName(obj,params)
% dirname = paramsToColorModulationDirName(obj,params)
%
% Generate a directory names that captures the background
% parameters

if (~strcmp(params.type,'Background'))
    error('Incorrect parameter type passed');
end

switch (params.backgroundType)
    case 'monitor'
        dirname = sprintf('lum%0.2f_lf%0.1f_leak%0.1f',...
            params.backgroundxyY(3),...
            params.lumFactor,...
            params.leakageLum);
    case 'AO'
        for ii = 1:length(params.backgroundWavelengthsNm)
            dirnameTmp = sprintf('%d_%0.1f',params.backgroundWavelengthsNm,params.backgroundCornealIrradianceUW);
            if (ii == 1)
                dirname = dirnameTmp;
            else
                dirname = [dirname '_' dirnameTmp];
            end
        end
    otherwise
        error('Unknown background type specified');
end





