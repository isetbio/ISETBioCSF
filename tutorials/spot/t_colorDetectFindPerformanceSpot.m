function validationData = t_colorDetectFindPerformanceSpot(varargin)
% validationData = t_colorDetectFindPerformanceSpot(varargin)
%
% This is a call into t_colorDetectFindPerformance that demonstates its
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
    
    % Override some defult parameters
    %
    % Set duration equal to sampling interval to do just one frame.
    rParams.temporalParams.stimulusDurationInSeconds = 200/1000;
    rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.stimulusDurationInSeconds;
    rParams.temporalParams.secondsToInclude = rParams.temporalParams.stimulusDurationInSeconds;

    % No eye movements
    rParams.temporalParams.emPathType = 'none';
    
    rParams.mosaicParams.integrationTimeInSeconds = rParams.temporalParams.stimulusDurationInSeconds;
    rParams.mosaicParams.isomerizationNoise = 'random';         % Type coneMosaic.validNoiseFlags to get valid values
    rParams.mosaicParams.osNoise = 'random';                    % Type outerSegment.validNoiseFlags to get valid values
    rParams.mosaicParams.osModel = 'Linear';
    
    rParams.oiParams.pupilDiamMm = 7;
end

%% Call into t_colorDetectFindPerformance with spot parameters
contrastParams = instanceParamsGenerate('instanceType','contrasts');
thresholdParams = thresholdParamsGenerate;
thresholdParams.method = 'mlpt';
validationData = t_colorDetectFindPerformance('rParams',rParams,'testDirectionParams',contrastParams,'thresholdParams',thresholdParams,'plotPsychometric',true,'generatePlots',p.Results.generatePlots);


