function dirname = paramsToResponseGenerationDirName(obj,oiParams)
% dirname = paramsToResponseGenerationDirName(obj,oiParams)
% 
% Generate a directory names that captures the oi parameters.

if (~strcmp(oiParams.type,'Optics'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('b%0.0f_l%0.0f_pup%0.1f', ...
    oiParams.blur, ...
    oiParams.lens, ...
    oiParams.pupilDiamMm);

