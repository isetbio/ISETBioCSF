function dirname = paramsToOiDirName(obj,oiParams)
% dirname = paramsToResponseGenerationDirName(obj,oiParams)
% 
% Generate a directory names that captures the oi parameters.

if (~strcmp(oiParams.type,'Optics'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('O_blur%0.0f_lens%0.0f_pupilDiam%0.1f_%s', ...
    oiParams.blur, ...
    oiParams.lens, ...
    oiParams.pupilDiamMm, ...
    oiParams.opticsModel); 

