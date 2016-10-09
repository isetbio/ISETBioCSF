function validationData = t_coneIsomerzationsMovieSpot(rParams)
% validationData = t_coneIsomerzationsMovieSpot(rParams)
%
% This is a call into t_coneIsomerizationsMovie that demonstates its
% ability to handle AO spots as well as Gabor modulations on monitors.

%% Clear
if (nargin == 0)
    ieInit; close all;
end

%% Fix random number generator so we can validate output exactly
rng(1);

%% Get the parameters we need
%
% responseParamsGenerate returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
if (nargin < 1 | isempty(rParams))
    rParams = responseParamsGenerate('spatialType','spot','backgroundType','AO','modulationType','AO');
    
    % Override some defaults to make more sense for our spot application
    rParams.oiParams.pupilDiamMm = 7;
end

%% Call into t_coneIsomerizationsMovie with spot parameters
validationData = t_coneIsomerizationsMovie(rParams);

%% Send back some validation data if requested
if (nargout > 0)
    validationData.maxIsomerizations = maxIsomerizations;
    validationData.minIsomerizations = minIsomerizations;
    validationData.contrasts = contrasts;
end

