function validationData = t_colorGaborConeCurrentEyeMovementsResponseInstancesOneFrame(rParams)
% validationData = t_colorGaborConeCurrentEyeMovementsResponseInstancesOneFrame(rParams)
%
% Show how to generate a number of response instances for a given stimulus
% condition, for a single stimulus frame without eye movements.
% This tutorial relies on routine
%   colorDetectResponseInstanceFastArrayConstruct
% which does most of the hard work.  The basic principles underlying colorDetectResponseInstanceFastArrayConstruct
% itself is demonstrated in tutorial 
%   t_colorGaborConeCurrentEyeMovementsMovie
% but the actual routine has some tricks to make it go fast.  There is also
% are routine
%   colorDetectResponseInstanceArrayConstruct
% that works more like the tutorial but is slower.
%
% This tutorial saves its output in a .mat file, which cah then read in by
%   t_colorGaborDetectFindPerformance
% which shows how to use the data to find the thresholds.
%
% The output goes into a place determined by
%   colorGaborDetectOutputDir
% which itself checks for a preference set by
%   ISETColorDetectPreferencesTemplate
% which you may want to edit before running this and other scripts that
% produce substantial output.  The output within the main output directory
% is sorted by directories whose names are computed from parameters.  This
% naming is done in routine
%   paramsToDirName.
%
% The returned validation structure allows this routine to be called from a
% validation script driven by the UnitTest toolbox.
%
% The tutorial produces output according to a scheme controlled by the
% specified IBIOColorDetect rwObject.
%
% See also:
%   t_colorGaborRespnseGenerationParams
%   t_colorGaborConeCurrentEyeMovementsMovie
%	t_colorGaborDetectFindPerformance
%   colorDetectResponseInstanceArrayConstruct
%   colorDetectResponseInstanceFastArrayConstruct

%% Clear
ieInit; clear; close all;

%% Fix random number generator so we can validate output exactly
rng(1);

%% Get the parameters we need
%
% t_colorGaborResponseGenerationParams returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
if (nargin < 1 | isempty(rParams))
    rParams = t_colorGaborResponseGenerationParams;
end

% Override some parameters
rParams.temporalParams.simulationTimeStepSecs = 200/1000;
rParams.temporalParams.stimulusDurationInSeconds = 0;
rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.simulationTimeStepSecs;
rParams.temporalParams.secondsToInclude = rParams.temporalParams.simulationTimeStepSecs;
rParams.temporalParams.eyesDoNotMove = true; 

rParams.mosaicParams.timeStepInSeconds = rParams.temporalParams.simulationTimeStepSecs;
rParams.mosaicParams.integrationTimeInSeconds = rParams.mosaicParams.timeStepInSeconds;
rParams.mosaicParams.isomerizationNoise = true;
rParams.mosaicParams.osNoise = true;
rParams.mosaicParams.osModel = 'Linear';

%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
rwObject.readProgram = '';
rwObject.writeProgram = mfilename;
rwObject.parentParamsList = {};

% These may only work on some computers, depending on what
% infrastructure is installed.
visualizeResponses = false;
exportToPDF = false;
renderVideo = false;

%% Parameters that define the LM instances we'll generate here
%
% Make these numbers small (trialNum = 2, deltaAngle = 180,
% nContrastsPerDirection = 2) to run through a test quickly.
LMPlaneInstanceParams = LMPlaneInstanceParamsGenerate;

%% Create the optics
theOI = colorDetectOpticalImageConstruct(rParams.oiParams);

%% Create the cone mosaic
theMosaic = colorDetectConeMosaicConstruct(rParams.mosaicParams);

%% Define stimulus set
%
% Chromatic directions in L/M plane.  It's a little easier to think in
% terms of angles.
LMangles = (0:LMPlaneInstanceParams.deltaAngle:180-LMPlaneInstanceParams.deltaAngle)/180*pi;
for angleIndex = 1:numel(LMangles)
    theta = LMangles(angleIndex);
    baseTestConeContrastDirs(:,angleIndex) = LMPlaneInstanceParams.baseStimulusLength*[cos(theta) sin(theta) 0.0]';
end

% Contrasts
if (strcmp(LMPlaneInstanceParams.contrastScale, 'linear'))
    testContrasts = linspace(LMPlaneInstanceParams.lowContrast, LMPlaneInstanceParams.highContrast, LMPlaneInstanceParams.nContrastsPerDirection);
else
    testContrasts = logspace(log10(LMPlaneInstanceParams.lowContrast), log10(LMPlaneInstanceParams.highContrast), LMPlaneInstanceParams.nContrastsPerDirection);
end

%% Generate data for the no stimulus condition
tic
colorModulationParamsTemp = rParams.colorModulationParams;
colorModulationParamsTemp.coneContrasts = [0 0 0]';
colorModulationParamsTemp.contrast = 0;
stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.2f', ...
    colorModulationParamsTemp.coneContrasts(1), colorModulationParamsTemp.coneContrasts(2), colorModulationParamsTemp.coneContrasts(3), colorModulationParamsTemp.contrast);
theNoStimData = struct(...
                 'testContrast', colorModulationParamsTemp.contrast, ...
            'testConeContrasts', colorModulationParamsTemp.coneContrasts, ...
                'stimulusLabel', stimulusLabel, ...
        'responseInstanceArray', colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, LMPlaneInstanceParams.trialsNum, rParams.temporalParams.simulationTimeStepSecs, ...
                                         rParams.gaborParams, colorModulationParamsTemp, rParams.temporalParams, theOI, theMosaic));
 
% Write the no cone contrast data
rwObject.currentParamsList = {rParams, colorModulationParamsTemp, LMPlaneInstanceParams};
rwObject.write('responseInstances_0',theNoStimData);

%% Generate data for all the examined stimuli
%
% This code is written in a slightly convoluted manner so that it will run
% inside a parfor loop -- we had to define some single indexed temp
% variables to get this to work.
%
% It is also possible that the parfor loop will not work for you, depending
% on your Matlab configuration.  In this case, change it to a for loop.

% Loop over color directions
parfor ii = 1:size(baseTestConeContrastDirs,2)
    % Find the highest in gamut cone contrast and define cone contrast
    % vector to be just under this length.
    gaborParamsLoop = rParams.gaborParams;
    gaborParamsLoop.coneContrasts = baseTestConeContrastDirs(:,ii);
    gaborParamsLoop.contrast = 1;
    [~,contrastScaleFactor(ii)] = colorGaborSceneCreate(gaborParamsLoop,true);
    testConeContrasts(:,ii) = 0.98*contrastScaleFactor(ii)*baseTestConeContrastDirs(:,ii);
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
        'responseInstanceArray', colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, LMPlaneInstanceParams.trialsNum, rParams.temporalParams.simulationTimeStepSecs, ...
                                          gaborParamsLoop, FOO, rParams.temporalParams, theOI, theMosaic));
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
                figHandle = visualizeResponseInstance(conditionDir, s.responseInstanceArray(iTrial), stimulusLabel, theMosaic, iTrial, LMPlaneInstanceParams.trialsNum, renderVideo);
                if (exportToPDF)
                    figFileNames{ii, jj, iTrial} = ...
                        fullfile(colorGaborDetectOutputDir(conditionDir),sprintf('%s_Trial%dOf%d.pdf', stimulusLabel, iTrial, LMPlaneInstanceParams.trialsNum),'figures');
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
