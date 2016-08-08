function dirname = paramsToResponseGenerationDirName(obj,params)
% pdirname = paramsToResponseGenerationDirName(obj,params)
% 
% Generate a directory names that captures the basic LMPlane instance generation
% parameters.

if (~strcmp(params.type,'LMPlaneInstance'))
    error('Incorrect parameter type passed');
end

theLMPlaneInstanceName = sprintf('LMPlane_ang%0.0f_ncon%0.0f',...
    params.deltaAngle,...
    params.nContrastsPerDirection ...
    );

dirname = theLMPlaneInstanceName;

