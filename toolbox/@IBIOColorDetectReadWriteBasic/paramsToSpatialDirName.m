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
        dirname = sprintf('ssz%0.0f_bsz%0.1f',...
            spotParams.spotSizeDegs,...
            spotParams.backgroundSizeDegs...
            );
    otherwise
        error('Unknown spatial type passed');
end