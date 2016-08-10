function dirname = paramsToGaborDirName(obj,gaborParams)
% dirname = paramsToGaborDirName(obj,gaborParams)
% 
% Generate a directory names that captures the basic non-color Gabor stimulus
% parameters.

if (~strcmp(gaborParams.type,'Gabor'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('cpd%0.0f_sfv%0.1f_fw%0.2f',...
    gaborParams.cyclesPerDegree,...
    gaborParams.fieldOfViewDegs,...
    gaborParams.gaussianFWHMDegs);