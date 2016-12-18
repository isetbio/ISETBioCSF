function validationData = t_coneIsomerizationsMovieSpot(varargin)
% validationData = t_coneIsomerizationsMovieSpot(varargin)
%
% This is a call into t_coneIsomerizationsMovie that demonstates its
% ability to handle AO spots as well as Gabor modulations on monitors.
%
% Optional key/value pairs
%  'rParams' - Value the is the rParams structure to use.  Default empty,
%     which then uses defaults produced by generation function.
%  'generatePlots' - true/fale (default true).  Make plots?

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('rParams',[],@isemptyorstruct);
p.addParameter('generatePlots',true,@islogical);
p.parse(varargin{:});
rParams = p.Results.rParams;

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
validationData = t_coneIsomerizationsMovie('rParams',rParams,'generatePlots',p.Results.generatePlots);


