function dirname = paramsToGaborDirName(obj,gaborParams)
% dirname = paramsToGaborDirName(obj,gaborParams)
% 
% Generate a directory names that captures the basic non-color stimulus
% parameters, as well as the oi and mosaic parameters used to generate the responses

if (~strcmp(gaborParams.type,'Gabor'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('cpd%0.0f_sfv%0.2f_fw%0.3f',...
    gaborParams.cyclesPerDegree,...
    gaborParams.fieldOfViewDegs,...
    gaborParams.gaussianFWHMDegs);