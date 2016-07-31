function dirname = paramsToResponseGenerationDirName(obj,params)
% dirname = paramsToResponseGenerationDirName(obj,params)
% 
% Generate a directory names that captures the basic LMPlane instance generation
% parameters.

if (~strcmp(params.type,'threshold'))
    error('Incorrect parameter type passed');
end

theThresholdName = sprintf('%s_%s_%dInt_kFold%d_PCA%d',...
    params.method, ...
    params.signalSource, ...
    params.nIntervals, ...
    params.kFold, ...
    params.PCAComponents ...
    );

dirname = theThresholdName;

