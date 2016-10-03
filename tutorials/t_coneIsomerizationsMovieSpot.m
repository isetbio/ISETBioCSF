function validationData = t_coneIsomerzationsMovieSpot(rParams)
% validationData = t_coneIsomerzationsMovieSpot(rParams)
%
% Illustrates the basic steps required to calculate cone isomerizations for
% a monochromatic spot on a background, where the key parameters that will
% vary are the size of the spot, the size of the background, and the
% radiometric properties of the spot and the background.  The reason we
% want to do this is so that we can make predictions for how thresholds
% will vary as we change these properties, particularly the size of the
% spot, for various radiometric, mosaic, and eccentricity choices.
%
% If parameters structure is not passed, the routine will use the defaults
% provided by
%   responseParamsGenerate('spatialType','spot','backgroundType','AO','modulationType','AO')
% That function and its subfunctions also documents what the relavant parameters are.
%
% The code illustrated here is encapsulated into function
%   colorSceneCreate.
%
% The returned validation structure allows this routine to be called from a
% validation script driven by the UnitTest toolbox.
%
% The tutorial produces output according to a scheme controlled by the
% specified IBIOColorDetect rwObject.
%
% See also:
%
% 9/14/16  dhb, wst  Wrote it.

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
end

%% Override some defaults to make more sense for our spot application
rParams.oiParams.pupilDiamMm = 7;

%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
theProgram = mfilename;
paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, rParams.colorModulationParams};

%% Call into t_coneIsomerizationsMovie with spot parameters
validationData = t_coneIsomerizationsMovie(rParams);

%% Send back some validation data if requested
if (nargout > 0)
    validationData.maxIsomerizations = maxIsomerizations;
    validationData.minIsomerizations = minIsomerizations;
    validationData.contrasts = contrasts;
end

