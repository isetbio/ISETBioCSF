function dirname = paramsToPsychoEllipsoidDirName(obj,params)
% dirname = paramsToPsychoEllipsoidDirName(obj,params)
% 
% Generate a directory names that captures the psychophysical parameters
% used to generate data/plots.

if (~strcmp(params.type,'psychoEllipsoid'))
    error('Incorrect parameter type passed');
end

theThresholdName = sprintf('%s_%s_sf%d',...
    params.subject, ...
    params.condition, ...
    params.theSf ...
    );

dirname = theThresholdName;

