function dirname = paramsToResponseGenerationDirName(obj,oiParams)
% dirname = paramsToResponseGenerationDirName(obj,oiParams)
% 
% Generate a directory names that captures the oi parameters.

if (~strcmp(oiParams.type,'Optics')) && (~strcmp(oiParams.type,'Optics_v2'))
    error('Incorrect parameter type passed');
end

if (strcmp(oiParams.type,'Optics'))
dirname = sprintf('b%0.0f_l%0.0f_pup%0.1f', ...
    oiParams.blur, ...
    oiParams.lens, ...
    oiParams.pupilDiamMm);
else
dirname = sprintf('[OPTICS]_blur%0.0f_lens%0.0f_pupilDiam%0.1f', ...
    oiParams.blur, ...
    oiParams.lens, ...
    oiParams.pupilDiamMm); 
end

