function validationData = t_colorDetectFindPerformanceSpot(varargin)
% Demonstrate t_colorDetectFindPerformance's ability to handle AO spots
%
% Syntax:
%   validationData = t_colorDetectFindPerformanceSpot(varargin)
%
% Description:
%    This is a call into t_colorDetectFindPerformance that demonstates its
%    ability to handle AO spots as well as Gabor modulations on monitors.
%
% Inputs:
%    None required.
%
% Outputs:
%    validationData - Struct. A structure of validation data.
%
% Optional key/value pairs:
%    rParams        - Struct. The rParams structure to use. Default empty,
%                     which uses defaults produced by generation function.
%    generatePlots  - Boolean. Whether to make plots. Default true.
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

    % Override some defult parameters
    %
    % Set duration equal to sampling interval to do just one frame.
    rParams.temporalParams.stimulusDurationInSeconds = 200 / 1000;
    rParams.temporalParams.stimulusSamplingIntervalInSeconds = ...
        rParams.temporalParams.stimulusDurationInSeconds;
    rParams.temporalParams.secondsToInclude = ...
        rParams.temporalParams.stimulusDurationInSeconds;

    % No eye movements
    rParams.temporalParams.emPathType = 'none';

    rParams.mosaicParams.integrationTimeInSeconds = ...
        rParams.temporalParams.stimulusDurationInSeconds;
    % For iso noise, type coneMosaic.validNoiseFlags to get valid values
    rParams.mosaicParams.isomerizationNoise = 'random';
    % For osNoise, type outerSegment.validNoiseFlags to get valid values
    rParams.mosaicParams.osNoise = 'random';
    rParams.mosaicParams.osModel = 'Linear';

    rParams.oiParams.pupilDiamMm = 7;
end

%% Call into t_colorDetectFindPerformance with spot parameters
contrastParams = instanceParamsGenerate('instanceType', 'contrasts');
thresholdParams = thresholdParamsGenerate;
thresholdParams.method = 'mlpt';
validationData = t_colorDetectFindPerformance('rParams', rParams, ...
    'testDirectionParams', contrastParams, ...
    'thresholdParams', thresholdParams, 'plotPsychometric', true, ...
    'generatePlots', p.Results.generatePlots);

end
