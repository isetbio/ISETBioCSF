function dirname = paramsToGaborDirName(obj,spatialParams)
% dirname = paramsToGaborDirName(obj,spatialParams)
% 
% Generate a directory names that captures the basic non-color Gabor stimulus
% parameters.

if (~strcmp(spatialParams.type,'Gabor'))
    error('Incorrect parameter type passed');
end

switch spatialParams.windowType
    case 'Gaussian';
        windowStr = sprintf('gfw%0.2f',spatialParams.gaussianFWHMDegs);
    case 'halfcos'
        windowStr = sprintf('cfw%0.2f',spatialParams.gaussianFWHMDegs);
end

dirname = sprintf('cpd%0.0f_sfv%0.1f_%s',...
    spatialParams.cyclesPerDegree,...
    spatialParams.fieldOfViewDegs,...
    windowStr);