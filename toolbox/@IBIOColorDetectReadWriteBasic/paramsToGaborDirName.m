function dirname = paramsToGaborDirName(obj,gaborParams)
% dirname = paramsToGaborDirName(obj,gaborParams)
% 
% Generate a directory names that captures the basic non-color Gabor stimulus
% parameters.

if (~strcmp(gaborParams.type,'Gabor'))
    error('Incorrect parameter type passed');
end

switch gaborParams.windowType
    case 'Gaussian';
        windowStr = sprintf('gfw%0.2f',gaborParams.gaussianFWHMDegs);
    case 'halfcos'
        windowStr = sprintf('cfw%0.2f',gaborParams.gaussianFWHMDegs);
end

dirname = sprintf('cpd%0.0f_sfv%0.1f_%s',...
    gaborParams.cyclesPerDegree,...
    gaborParams.fieldOfViewDegs,...
    windowStr);