function validationData = t_coneIsomerizationsMovie(rParams,varargin)
% validationData = t_coneIsomerizationsMovie(rParams,varargin)
%
% Illustrates the basic steps required to calculate cone isomerizations
% for a Gaussian windowed temporal color modulation.
%
% If parameters structure is not passed, the routine will use the defaults
% provided by
%   responseParamsGenerate
% This will produce a Gabor pattern.
% That function's subfunctions also document what the relavant parameters are.
%
% t_coneIsomerizationsMovieSpot calls into this one with parameters that
% produce a demo for a spot.
%
% The returned validation structure allows this routine to be called from a
% validation script driven by the UnitTest toolbox.
%
% The tutorial produces output according to a scheme controlled by the
% specified IBIOColorDetect rwObject.
%
% See also:  
%   t_coneIsomerizationsMovieSpot
%	t_coneCurrentEyeMovementsMovie
%   responseParamsGenerate
%   colorSceneCreate
%
% Optional key/value pairs
%  'isomerizationNoise' - true/false (default false).  Add isomerization noise?
%  'generatePlots' - true/false (default true).  Make plots?
%  'setRngSeed' - true/false (default true). When true, set the rng seed so noise is frozen.

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('isomerizationNoise','none', @(x)ismember(x, coneMosaic.validNoiseFlags));
p.addParameter('generatePlots',true,@islogical);
p.addParameter('setRngSeed',true,@islogical);
p.parse(varargin{:});

%% Clear
if (nargin == 0)
    ieInit; close all;
end

%% Fix random number generator so we can validate output exactly
if (p.Results.setRngSeed)
    rng(1);
end;

%% Get the parameters we need
if (nargin < 1 | isempty(rParams))
    rParams = responseParamsGenerate;
    
    % Override some of the defaults
    rParams.mosaicParams.isomerizationNoise = p.Results.isomerizationNoise;
end

%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
theProgram = mfilename;
paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, rParams.colorModulationParams};

%% Plot the Gaussian temporal window, just to make sure it looks right
if (p.Results.generatePlots)
    gaussianFigure = figure; clf;
    plot(rParams.temporalParams.sampleTimes,rParams.temporalParams.gaussianTemporalWindow,'r');
    xlabel('Time (seconds)');
    ylabel('Window Amplitude');
    title('Stimulus Temporal Window');
end



%% Loop over time and build a cell array of scenes
gaborScene = cell(rParams.temporalParams.nSampleTimes,1);
for ii = 1:rParams.temporalParams.nSampleTimes
    % Make the scene for this timestep.  We make a temporary parameters
    % structure and tweak the contrast according to the temporal Gaussian
    % window.
    rParamsTemp = rParams;
    rParamsTemp.colorModulationParams.contrast = rParams.colorModulationParams.contrast*rParams.temporalParams.gaussianTemporalWindow(ii);
    fprintf('Computing scene %d of %d, time %0.3f, windowVal %0.3f\n',ii,rParamsTemp.temporalParams.nSampleTimes,rParamsTemp.temporalParams.sampleTimes(ii),rParamsTemp.temporalParams.gaussianTemporalWindow(ii));
    gaborScene{ii} = colorSceneCreate(rParamsTemp.spatialParams,rParams.backgroundParams,rParamsTemp.colorModulationParams,rParams.oiParams);
end
clearvars('rParamsTemp');

% Make a movie of the stimulus sequence
if (p.Results.generatePlots)
    showLuminanceMap = false;
    visualizeSceneOrOpticalImageSequence(rwObject,paramsList,theProgram, ...
        'scene', gaborScene, rParams.temporalParams.sampleTimes, showLuminanceMap, 'gaborStimulusMovie');
end

%% Create the OI object we'll use to compute the retinal images from the scenes
%
% Then loop over scenes and compute the optical image for each one 
theBaseOI = colorDetectOpticalImageConstruct(rParams.oiParams);
theOI = cell(rParams.temporalParams.nSampleTimes,1);
for ii = 1:rParams.temporalParams.nSampleTimes
    % Compute retinal image
    fprintf('Computing optical image %d of %d, time %0.3f\n',ii,rParams.temporalParams.nSampleTimes,rParams.temporalParams.sampleTimes(ii));
    theOI{ii} = oiCompute(theBaseOI,gaborScene{ii});
end

% Make a movie of the optical image sequence
if (p.Results.generatePlots)
showLuminanceMap = false;
visualizeSceneOrOpticalImageSequence(rwObject,paramsList,theProgram, ...
    'optical image', theOI, rParams.temporalParams.sampleTimes, showLuminanceMap, 'gaborOpticalImageMovie');
end

%% Create the coneMosaic object we'll use to compute cone respones
rParams.mosaicParams.fieldOfViewDegs = rParams.spatialParams.fieldOfViewDegs;
theMosaic = colorDetectConeMosaicConstruct(rParams.mosaicParams);
for ii = 1:rParams.temporalParams.nSampleTimes      
    % Compute mosaic response for each stimulus frame
    % For a real calculation, we would save these so that we could use them
    % to do something.  But here we just (see next line) compute the
    % contrast seen by each class of cone in the mosaic, just to show we
    % can do something.
    fprintf('Computing isomerizations %d of %d, time %0.3f\n',ii,rParams.temporalParams.nSampleTimes,rParams.temporalParams.sampleTimes(ii));
    gaborConeAbsorptions(:,:,ii) = theMosaic.compute(theOI{ii},'currentFlag',false);
end

% Make a movie of the isomerizations
if (p.Results.generatePlots)
    eyeMovementSequence = [];
    visualizeMosaicResponseSequence(rwObject,paramsList,theProgram, ...
        'isomerizations (R*/cone)', gaborConeAbsorptions, eyeMovementSequence, theMosaic.pattern, rParams.temporalParams.sampleTimes, [theMosaic.width theMosaic.height], theMosaic.fov, rParams.mosaicParams.integrationTimeInSeconds, 'gaborIsomerizations');
end

%% Plot cone contrasts as a function of time, as a check
%
% As with the contrats in t_colorGabor, these are very close to right,
% although not completely perfect.
for ii = 1:rParams.temporalParams.nSampleTimes      
    LMSContrasts(:,ii) = mosaicUnsignedConeContrasts(gaborConeAbsorptions(:,:,ii),theMosaic);
end
if (p.Results.generatePlots)
    vcNewGraphWin; hold on;
    plot(rParams.temporalParams.sampleTimes,LMSContrasts(1,:)','r');
    plot(rParams.temporalParams.sampleTimes,LMSContrasts(2,:)','g');
    plot(rParams.temporalParams.sampleTimes,LMSContrasts(3,:)','b');
    xlabel('Time (seconds)');
    ylabel('Contrast');
    title('LMS Cone Contrasts');
end

%% Send back some validation data if requested
if (nargout > 0)
    validationData.LMSContrats = LMSContrasts;
end

