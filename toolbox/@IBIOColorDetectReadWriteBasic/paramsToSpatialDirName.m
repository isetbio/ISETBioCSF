function dirname = paramsToSpatialDirName(obj,spatialParams)
% dirname = paramsToSpatialDirName(obj,spatialParams)
%
% Generate a directory names that captures the basic spatial stimulus
% parameters.

if (~strcmp(spatialParams.type,'Spatial')) && (~strcmp(spatialParams.type,'Spatial_v2'))
    error('Incorrect parameter type passed');
end

if (strcmp(spatialParams.type,'Spatial'))
    switch (spatialParams.spatialType)
        case 'Gabor'
            switch spatialParams.windowType
                case 'Gaussian';
                    windowStr = sprintf('gfw%0.2f',spatialParams.gaussianFWHMDegs);
                case 'halfcos'
                    windowStr = sprintf('cfw%0.2f',spatialParams.gaussianFWHMDegs);
                otherwise
                    error('Unknown window type field passed');
            end
            dirname = sprintf('cpd%0.0f_sfv%0.1f_%s',...
                spatialParams.cyclesPerDegree,...
                spatialParams.fieldOfViewDegs,...
                windowStr);
        case 'spot'
            dirname = sprintf('ssz%0.3f_bsz%0.1f_sfv%0.1f',...
                spatialParams.spotSizeDegs,...
                spatialParams.backgroundSizeDegs,...
                spatialParams.fieldOfViewDegs...                 
                );
        otherwise
            error('Unknown spatial type passed');
    end
else
    switch (spatialParams.spatialType)
        case 'Gabor'
            switch spatialParams.windowType
                case 'Gaussian';
                    windowStr = sprintf('Gaussian_FWHMDegs%0.2f',spatialParams.gaussianFWHMDegs);
                case 'halfcos'
                    windowStr = sprintf('HalfCos_FWHMDegs%0.2f',spatialParams.gaussianFWHMDegs);
                otherwise
                    error('Unknown window type field passed');
            end
            dirname = sprintf('[STIM_SPATIAL]_SFcpd%0.2f_FOVdegs%0.1f_%s',...
                spatialParams.cyclesPerDegree,...
                spatialParams.fieldOfViewDegs,...
                windowStr);
        otherwise
            error('Unknown spatial type passed');
    end
            
end
