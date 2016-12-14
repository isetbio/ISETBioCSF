function [validationData, extraData] = t_coneCurrentEyeMovementsResponseInstancesSpot(varargin)
% [validationData, extraData] = t_coneCurrentEyeMovementsResponseInstancesSpot(varargin)
%
% This is a call into t_coneCurrentEyeMovementsResponseInstances that demonstates its
% ability to handle AO spots as well as Gabor modulations on monitors.
%
% Key/value pairs
%   'rParams' - structure (default empty). Value the is the rParams structure to use
%   'contrastParams - structure (default empty). Value is the contrastParams structure to use.
%   'freezeNoise' - true/false (default true).  Freezes all noise so that results are reproducible
%   'compute' - true/false (default true).  Do the computations.
%   'generatePlots' - true/false (default false).  Produce response
%        visualizations.  Set to false when running big jobs on clusters or
%        in parfor loops, as plotting doesn't seem to play well with those
%        conditions.
%   'exportPDF' - true/false (default true).  If visualizing responses,
%        export the PDF files.
%   'renderVideo' - true/false (default true).  If visualizing responses, generate
%        the videos.
%   'delete' - true/false (default true).  Delete the response instance
%        files.  Useful for cleaning up big output when we are done with
%        it.  If this is true, output files are deleted at the end.

%% Parse input
p = inputParser;
p.addParameter('rParams',[],@isemptyorstruct);
p.addParameter('contrastParams',[],@isemptyorstruct);
p.addParameter('freezeNoise',true,@islogical);
p.addParameter('compute',true,@islogical);
p.addParameter('generatePlots',false,@islogical);
p.addParameter('exportPDF',true,@islogical);
p.addParameter('renderVideo',true,@islogical);
p.addParameter('delete',false',@islogical);
p.parse(varargin{:});
rParams = p.Results.rParams;
contrastParams = p.Results.contrastParams;

%% Clear
if (nargin == 0)
    ieInit; close all;
end

%% Get the parameters we need
%
% responseParamsGenerate returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
if (isempty(rParams))
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

% Fix random number generator so we can validate output exactly
if (p.Results.freezeNoise)
     fprintf(2, '\n%s: freezing all noise \n', mfilename);
     rng(1);
     if (strcmp(rParams.mosaicParams.isomerizationNoise, 'random'))
         fprintf(2, '\tmosaicParams.isomerizationNoise was set to ''%s'', setting it to ''frozen''.\n', rParams.mosaicParams.isomerizationNoise);
         rParams.mosaicParams.isomerizationNoise = 'frozen';
     end
     if (strcmp(rParams.mosaicParams.osNoise, 'random'))
         fprintf(2, '\tmosaicParams.osNoise was set to ''%s'', setting it to ''frozen''.\n', rParams.mosaicParams.osNoise);
         rParams.mosaicParams.osNoise = 'frozen';
     end
end
 
%% Call into t_coneCurrentEyeMovementsResponseInstances with spot parameters
if (isempty(contrastParams))
    contrastParams = instanceParamsGenerate('instanceType','contrasts');
end
[validationData, extraData] = t_coneCurrentEyeMovementsResponseInstances('rParams',rParams,'testDirectionParams',contrastParams, 'freezeNoise', p.Results.freezeNoise, 'generatePlots',p.Results.generatePlots);


