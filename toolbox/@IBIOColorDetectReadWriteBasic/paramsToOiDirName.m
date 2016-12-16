function dirname = paramsToOiDirName(obj,oiParams)
% dirname = paramsToResponseGenerationDirName(obj,oiParams)
% 
% Generate a directory names that captures the oi parameters.

if (~strcmp(oiParams.type,'Optics'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('[OPTICS]_blur%0.0f_lens%0.0f_pupilDiam%0.1f', ...
    oiParams.blur, ...
    oiParams.lens, ...
    oiParams.pupilDiamMm); 

