%% t_colorGaborConeCurrentEyeMovementsResponseInstances
%
% Show how to generate a number of response instances for a given stimulus condition.
% This tutorial relies on routine
%   colorDetectResponseInstanceArrayConstruct
% which does most of the hard work.  The basic code underlying colorDetectResponseInstanceArrayConstruct
% itself is demonstrated in tutorial 
%   t_colorGaborConeCurrentEyeMovementsMovie.
%
% This tutorial saves its output in a .mat file, which is then read in by
%   t_colorGaborDetectFindThresholds
% which shows how to use the data to find the thresholds.
%
% 7/9/16  npc Wrote it.

%% Initialize
ieInit; clear; close all;

% Add project toolbox to Matlab path
AddToMatlabPathDynamically(fullfile(fileparts(which(mfilename)),'../toolbox')); 

%% Parameters that control output

% Set to true to save data for use by t_colorGaborDetectFindThresholds
saveData = true;

% These may only work on some computers, depending on what
% infrastructure is installed.
visualizeResponses = false;
exportToPDF = false;
renderVideo = false;

%% Parameters that define how much we do here

% Define how many noisy data instances to generate
trialsNum = 500; %500;

% Delta angle sampling in LM plane (samples between 0 and 180 degrees)
%
% Also base stimulus length in cone contrast space.  This variable
% no long has an effect because we scale each base direction to be 
% just inside monitor gamut
deltaAngle = 15; % 15; 
baseStimulusLength = 1;

% Number of contrasts to run in each color direction
nContrastsPerDirection = 10; % 10;
lowContrast = 0.001;
highContrast = 0.4;

%% Define parameters of simulation
%
% The time step at which to compute eyeMovements and osResponses
simulationTimeStep = 10/1000;

% Stimulus (gabor) params
gaborParams.fieldOfViewDegs = 1.0;
gaborParams.gaussianFWHMDegs = 0.35;
gaborParams.cyclesPerDegree = 2;
gaborParams.row = 128;
gaborParams.col = 128;
gaborParams.ang = 0;
gaborParams.ph = 0;
gaborParams.backgroundxyY = [0.27 0.30 49.8]';
gaborParams.leakageLum = 2.0;
gaborParams.monitorFile = 'CRT-MODEL';
gaborParams.viewingDistance = 0.75;

% Temporal modulation and stimulus sampling parameters.
%
% The millisecondsToInclude field tells how many milliseconds of the
% stimulus around the stimulus peak to include in data saved to pass to the
% classification routines.  Since the response might be delayed, we also
% allow specification of an offset (positive numbers mean later) around
% which the window is taken.
frameRate = 60;
temporalParams.windowTauInSeconds = 0.165;
temporalParams.stimulusDurationInSeconds = 2*temporalParams.windowTauInSeconds;
temporalParams.stimulusSamplingIntervalInSeconds = 1/frameRate;
temporalParams.millisecondsToInclude = 50;
temporalParams.millisecondsToIncludeOffset = 35;

% Optionally, have zero amplitude eye movements
temporalParams.eyesDoNotMove = false; 

% Optional CRT raster effects.
% 
% The underlying routine that generates temporal samples 
% can simulate the fact that CRTs produce an impulse during
% each frame, although this simulation works on a frame basis
% not on a pixel-by-pixel basis.  
% 
% The parameer rasterSamples is the number
% of raster samples generated per CRT refresh
% interval.
temporalParams.addCRTrasterEffect = false;
temporalParams.rasterSamples = 5; 
if (temporalParams.addCRTrasterEffect)
    simulationTimeStep = simulationTimeStep/temporalParams.rasterSamples;
end

% Optical image parameters
oiParams.fieldOfViewDegs = gaborParams.fieldOfViewDegs;
oiParams.offAxis = false;
oiParams.blur = false;
oiParams.lens = true;

% Cone mosaic parameters
mosaicParams.fieldOfViewDegs = gaborParams.fieldOfViewDegs;
mosaicParams.macular = true;
mosaicParams.LMSRatio = [0.62 0.31 0.07];
mosaicParams.timeStepInSeconds = simulationTimeStep;
mosaicParams.integrationTimeInSeconds = mosaicParams.timeStepInSeconds;
mosaicParams.photonNoise = true;
mosaicParams.osNoise = true;
mosaicParams.osModel = 'Linear';

%% Create the optics
theOI = colorDetectOpticalImageConstruct(oiParams);

%% Create the cone mosaic
theMosaic = colorDetectConeMosaicConstruct(mosaicParams);

%% Define stimulus set
%
% Chromatic directions in L/M plane.  It's a little easier to think in
% terms of angles.
LMangles = (0:deltaAngle:180-deltaAngle)/180*pi;
for angleIndex = 1:numel(LMangles)
    theta = LMangles(angleIndex);
    testConeContrastsDirs(:,angleIndex) = baseStimulusLength*[cos(theta) sin(theta) 0.0]';
end

% Contrasts
testContrasts = linspace(lowContrast, highContrast, nContrastsPerDirection);

%% Generate data for the no stimulus condition
tic
gaborParams.coneContrasts = [0 0 0]';
gaborParams.contrast = 0;
stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.2f', gaborParams.coneContrasts(1), gaborParams.coneContrasts(2), gaborParams.coneContrasts(3), gaborParams.contrast);
theNoStimData = struct(...
                 'testContrast', gaborParams.contrast, ...
            'testConeContrasts', gaborParams.coneContrasts, ...
                'stimulusLabel', stimulusLabel, ...
        'responseInstanceArray', colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, trialsNum, simulationTimeStep, ...
                                         gaborParams, temporalParams, theOI, theMosaic));
 
%% Set up output dir
conditionDir = paramsToDirName(gaborParams,temporalParams,oiParams,mosaicParams,[]);
outputDir = colorGaborDetectOutputDir(conditionDir,'output');

%% Generate data for all the examined stimuli
%
% This code is written in a slightly convoluted manner so that it will run
% inside a parfor loop -- we had to define some single indexed temp
% variables to get this to work.
%
% It is also possible that the parfor loop will not work for you, depending
% on your Matlab configuration.  In this case, change it to a for loop.

% Loop over color directions
parfor ii = 1:size(testConeContrastsDirs,2)
    % Find the highest in gamut cone contrast and define cone contrast
    % vector to be just under this length.
    gaborParamsLoop = gaborParams;
    gaborParamsLoop.coneContrasts = testConeContrastsDirs(:,ii);
    gaborParamsLoop.contrast = 1;
    [~,contrastScaleFactor(ii)] = colorGaborSceneCreate(gaborParamsLoop,true);
    testConeContrasts(:,ii) = 0.98*contrastScaleFactor(ii)*testConeContrastsDirs(:,ii);
    gaborParamsLoop.coneContrasts = testConeContrasts(:,ii);
    
    % Make noisy instances for each contrast
    stimDataJJ = cell(1,numel(testContrasts));
    for jj = 1:numel(testContrasts)
        gaborParamsLoop.contrast = testContrasts(jj);
        stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.2f',...
            gaborParamsLoop.coneContrasts(1), gaborParamsLoop.coneContrasts(2), gaborParamsLoop.coneContrasts(3), gaborParamsLoop.contrast);
        stimDataJJ{jj} = struct(...
                 'testContrast', gaborParamsLoop.contrast, ...
            'testConeContrasts', gaborParamsLoop.coneContrasts, ...
                'stimulusLabel', stimulusLabel, ...
        'responseInstanceArray', colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, trialsNum, simulationTimeStep, ...
                                          gaborParamsLoop, temporalParams, theOI, theMosaic));
    end
    
    % Save data for this color direction
    parforResponseInstancesSave(stimDataJJ,fullfile(outputDir,sprintf('responseInstances_%d',ii)));
end 
fprintf('Finished generating responses in %2.2f minutes\n', toc/60);

%% Save the other data we need for use by the classifier preprocessing subroutine
%
% And also a copy of this script
save(fullfile(outputDir,'responseInstances_0'), 'theNoStimData', 'testConeContrasts', 'testContrasts', 'theMosaic', 'gaborParams', 'temporalParams', 'oiParams', 'mosaicParams', '-v7.3');
scriptDir = colorGaborDetectOutputDir(conditionDir,'scripts');
unix(['cp ' mfilename('fullpath') '.m ' scriptDir]);

%% Visualize responses
%
% PROBABLY BROKEN: This will probably not work now that we have mucked with
% the way the data get saved.
%
% Also, the time numbers on the videos do not seem to correspond to the
% stimulus peak at time 0, and the videos look screwy in the no noise
% condition.
if (visualizeResponses)
    fprintf('\nVisualizing responses ...\n');
    for ii = 1:size(testConeContrasts,2)
        for jj = 1:numel(testContrasts)
            stimulusLabel = sprintf('LMS_%2.2f_%2.2f_%2.2f_Contrast_%2.2f', testConeContrasts(1,ii), testConeContrasts(2,ii), testConeContrasts(3,ii), testContrasts(jj));
            s = theStimData{ii, jj}; 
            
            % Visualize a few response instances only
            for iTrial = 1:2
                figHandle = visualizeResponseInstance(conditionDir, s.responseInstanceArray(iTrial), stimulusLabel, theMosaic, iTrial, trialsNum, renderVideo);
                if (exportToPDF)
                    figFileNames{ii, jj, iTrial} = ...
                        fullfile(colorGaborDetectOutputDir(conditionDir),sprintf('%s_Trial%dOf%d.pdf', stimulusLabel, iTrial, trialsNum),'figures');
                    NicePlot.exportFigToPDF(figFileNames{ii, jj, iTrial}, figHandle, 300);
                end
            end 
        end
    end

    % Export summary PDF with all responses
    if (exportToPDF)
        summaryPDF = fullfile(colorGaborDetectOutputDir(conditionDir,'figures'), 'AllInstances.pdf');
        fprintf('Exporting a summary PDF with all response instances in %s\n', summaryPDF);
        NicePlot.combinePDFfilesInSinglePDF(figFileNames(:), summaryPDF);
    end
end
