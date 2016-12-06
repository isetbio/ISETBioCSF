function validationData = t_coneCurrentEyeMovementsResponseInstancesSpot(varargin)
% validationData = t_coneCurrentEyeMovementsResponseInstancesSpot(varargin)
%
% This is a call into t_coneCurrentEyeMovementsResponseInstances that demonstates its
% ability to handle AO spots as well as Gabor modulations on monitors.
%
% Key/value pairs
%   'rParams' - structure (default empty). Value the is the rParams structure to use
%   'contrastParams - structure (default empty). Value is the contrastParams structure to use.
%   'setRngSeed' - true/false (default true).  Set the rng seed to a
%        value so output is reproducible.
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
p.addParameter('setRng',true,@islogical);
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

%% Fix random number generator so we can validate output exactly
rng(1);

%% Get the parameters we need
%
% responseParamsGenerate returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
if (isempty(rParams))
    rParams = responseParamsGenerate('spatialType','spot','backgroundType','AO','modulationType','AO');
    
    % Override some defult parameters
    %
    % Set duration equal to sampling interval to do just one frame.
    rParams.temporalParams.simulationTimeStepSecs = 200/1000;
    rParams.temporalParams.stimulusDurationInSeconds = rParams.temporalParams.simulationTimeStepSecs;
    rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.simulationTimeStepSecs;
    rParams.temporalParams.secondsToInclude = rParams.temporalParams.simulationTimeStepSecs;
    rParams.temporalParams.eyesDoNotMove = true;
    
    rParams.mosaicParams.timeStepInSeconds = rParams.temporalParams.simulationTimeStepSecs;
    rParams.mosaicParams.integrationTimeInSeconds = rParams.mosaicParams.timeStepInSeconds;
    rParams.mosaicParams.isomerizationNoise = 'random';        % select from {'random', 'frozen', 'none'}
    rParams.mosaicParams.osNoise = 'random';        % select from {'random', 'frozen', 'none'}
    rParams.mosaicParams.osModel = 'Linear';
    
    rParams.oiParams.pupilDiamMm = 7;
end

%% Call into t_coneCurrentEyeMovementsResponseInstances with spot parameters
if (isempty(contrastParams))
    contrastParams = instanceParamsGenerate('instanceType','contrasts');
end
validationData = t_coneCurrentEyeMovementsResponseInstances('rParams',rParams,'testDirectionParams',contrastParams,'generatePlots',p.Results.generatePlots);


