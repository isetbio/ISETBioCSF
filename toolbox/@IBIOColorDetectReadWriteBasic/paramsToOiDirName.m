function dirname = paramsToResponseGenerationDirName(obj,oiParams)
% dirname = paramsToResponseGenerationDirName(obj,oiParams)
% 
% Generate a directory names that captures the basic non-color stimulus
% parameters, as well as the oi and mosaic parameters used to generate the responses

if (~strcmp(oiParams.type,'Optics'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('b%0.0f_l%0.0f', ...
    oiParams.blur, ...
    oiParams.lens);

