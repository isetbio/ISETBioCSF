function dirname = paramsToResponseGenerationDirName(obj,params)
% pdirname = paramsToResponseGenerationDirName(obj,params)
% 
% Generate a directory names that captures the basic LMPlane instance generation
% parameters.

if (~strcmp(params.type,'LMPlaneInstance'))
    error('Incorrect parameter type passed');
end

theLMPlaneInstanceName = sprintf('LMPlane_sang%0.0f_dang%0.0f_nang%0.0f_ncon%0.0f_nt%0.0f',...
    params.startAngle,...
    params.deltaAngle,...
    params.nAngles,...
    params.nContrastsPerDirection, ...
    params.trialsNum ...
    );

dirname = theLMPlaneInstanceName;

