function validationData = t_coneIsomerizationsMovieSpot(varargin)
% Demonstrate how to handle AO spots and Gabor modulations on monitors.
%
% Syntax:
%   validationData = t_coneIsomerizationsMovieSpot([varargin])
%
% Description:
%    This is a call into t_coneIsomerizationsMovie that demonstates its
%    ability to handle AO spots as well as Gabor modulations on monitors.
%
% Inputs:
%    None required.
%
% Outputs:
%    validationData - Struct. A structure containing the validation data.
%
% Optional key/value pairs
%    rParams        - Struct. Structure of parameters. Default []. Default
%                     will then call the generation function.
%    generatePlots  - Boolean. Whether to generate plots. Default true.
%

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('rParams', [], @isemptyorstruct);
p.addParameter('generatePlots', true, @islogical);
p.parse(varargin{:});
rParams = p.Results.rParams;

%% Clear
if (nargin == 0)
    ieInit;
    close all;
end

%% Fix random number generator so we can validate output exactly
rng(1);

%% Get the parameters we need
% responseParamsGenerate returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
if (nargin < 1 | isempty(rParams))
    rParams = responseParamsGenerate('spatialType', 'spot', ...
        'backgroundType', 'AO', 'modulationType', 'AO');

    % Override some defaults to make more sense for our spot application
    rParams.oiParams.pupilDiamMm = 7;
end

%% Call into t_coneIsomerizationsMovie with spot parameters
validationData = t_coneIsomerizationsMovie('rParams', rParams, ...
    'generatePlots', p.Results.generatePlots);

end
