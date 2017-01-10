function dirname = paramsToSpatialDirName(obj,spatialParams)
% dirname = paramsToSpatialDirName(obj,spatialParams)
%
% Generate a directory names that captures the basic spatial stimulus
% parameters.

if (~strcmp(spatialParams.type,'Spatial'))
    error('Incorrect parameter type passed');
end

switch (spatialParams.spatialType)
    case 'Gabor'
        switch spatialParams.windowType
            case 'Gaussian'
                windowStr = sprintf('Gaussian_FWHMDegs%0.2f',spatialParams.gaussianFWHMDegs);
            case 'halfcos'
                windowStr = sprintf('HalfCos_FWHMDegs%0.2f',spatialParams.gaussianFWHMDegs);
            otherwise
                error('Unknown window type field passed');
        end
        dirname = sprintf('[STIM_SPATIAL]_SFcpd%0.2f_phase%0.2f_FOVdegs%0.2f_%s',...
            spatialParams.cyclesPerDegree,...
            spatialParams.ph,....
            spatialParams.fieldOfViewDegs,...
            windowStr);
    case 'spot'
        dirname = sprintf('[STIM_SPATIAL]_spotSideDegs%0.3f_backgroundSizeDegs%0.2f_fielfOfViewDegs%0.2f',...
            spatialParams.spotSizeDegs,...
            spatialParams.backgroundSizeDegs,...
            spatialParams.fieldOfViewDegs...                 
            );
    otherwise
        error('Unknown spatial type passed');
end
