function [validationData, extraData] = ...
    t_coneCurrentEyeMovementsResponseInstancesSpot(varargin)
% Demonstrate AO spots using t_coneCurrentEyeMovementsResponseInstances
%
% Syntax:
%   [validationData, extraData] = ...
%       t_coneCurrentEyeMovementsResponseInstancesSpot([varargin])
%
% Description:
%    This is a call into t_coneCurrentEyeMovementsResponseInstances that
%    demonstates its ability to handle AO spots as well as Gabor
%    modulations on monitors.
%
% Inputs:
%    None required.
%
% Outputs:
%    validationData - Struct. A structure of validation data.
%    extraData      - Struct. A structure of additional relevant data.
%
% Optional key/value pairs:
%    employStandardHostComputerResources
%                   - Boolean. Only validation scripts should set this to
%                     true, so as to produce identical sequences of random
%                     numbers. All other scripts should set it (leave it)
%                     to false. Default false.
%    rParams        - Struct. The rParams structure to use. Default empty,
%                     which uses defaults produced by generation function.
%    contrastParams - Struct. A structure of the contrastParams to use.
%                     Default empty, which takes advantage of a generation
%                     function to retrieve values.
%    freezeNoise    - Boolean. Whether to freeze all noise so that results
%                     are reproducible. Default true.
%    compute        - Boolean. Whether or not to perform the computations.
%                     Default true.
%    generatePlots  - Boolean. Whether to produce response visualizations.
%                     Set to false when running big jobs on clusters or in
%                     parfor loops, as plotting doesn't seem to play well
%                     with those conditions. Default false.
%    exportPDF      - Boolean. If visualizaing responses, export the PDF
%                     files. Default true.
%    renderVideo    - Boolean. If visualizing responses, generate the
%                     corresponding videos. Default true.
%    delete         - Boolean. Whether to delete the respones instance
%                     files. This is useful for cleaning up big output when
%                     we are done with it. If this is true, output files
%                     are deleted at the end. Default true.

%% Parse input
p = inputParser;
p.addParameter('rParams', [], @isemptyorstruct);
p.addParameter('employStandardHostComputerResources', false, @islogical);
p.addParameter('contrastParams', [], @isemptyorstruct);
p.addParameter('freezeNoise', true, @islogical);
p.addParameter('compute', true, @islogical);
p.addParameter('generatePlots', false, @islogical);
p.addParameter('exportPDF', true, @islogical);
p.addParameter('renderVideo', true, @islogical);
p.addParameter('delete', false', @islogical);
p.parse(varargin{:});
rParams = p.Results.rParams;
contrastParams = p.Results.contrastParams;

%% Clear
if (nargin == 0)
    ieInit;
    close all;
end

%% Get the parameters we need
% responseParamsGenerate returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
if (isempty(rParams))
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

% Fix random number generator so we can validate output exactly
if (p.Results.freezeNoise)
     fprintf(2, '\n%s: freezing all noise \n', mfilename);
     rng(1);
     if (strcmp(rParams.mosaicParams.isomerizationNoise, 'random'))
         fprintf(2, ['\tmosaicParams.isomerizationNoise was set to ', ...
             '''%s'', setting it to ''frozen''.\n'], ...
             rParams.mosaicParams.isomerizationNoise);
         rParams.mosaicParams.isomerizationNoise = 'frozen';
     end
     if (strcmp(rParams.mosaicParams.osNoise, 'random'))
         fprintf(2, ['\tmosaicParams.osNoise was set to ''%s'', ', ...
             'setting it to ''frozen''.\n'], rParams.mosaicParams.osNoise);
         rParams.mosaicParams.osNoise = 'frozen';
     end
end

%% Call t_coneCurrentEyeMovementsResponseInstances with spot parameters
if (isempty(contrastParams))
    contrastParams = instanceParamsGenerate('instanceType', 'contrasts');
end

[validationData, extraData] = ...
    t_coneCurrentEyeMovementsResponseInstances(...
    'employStandardHostComputerResources', ...
    p.Results.employStandardHostComputerResources, ...
    'rParams', rParams, 'testDirectionParams', contrastParams, ...
    'freezeNoise', p.Results.freezeNoise, ...
    'generatePlots', p.Results.generatePlots);

end
