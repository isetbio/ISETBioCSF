function validationData = t_colorGaborConeCurrentEyeMovementsResponseInstancesOneFrame(rParams,LMPlaneInstanceParams)
% validationData = t_colorGaborConeCurrentEyeMovementsResponseInstancesOneFrame([rParams],[LMPlaneInstanceParams])
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
    rParams = colorGaborResponseParamsGenerate;
    
    % Override some defult parameters
    %
    % A stimulus duration of 0 denotes a single frame static output
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
end

%% Parameters that define the LM instances we'll generate here
%
% Make these numbers small (trialNum = 2, deltaAngle = 180,
% nContrastsPerDirection = 2) to run through a test quickly.
if (nargin < 2 | isempty(LMPlaneInstanceParams))
    LMPlaneInstanceParams = LMPlaneInstanceParamsGenerate;
end

%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
theProgram = mfilename;
parentParamsList = {};
currentParamsList = {rParams, rParams.colorModulationParams};

%% Control what types of output are written
% These may only work on some computers, depending on what infrastructure is installed.
visualizeResponses = false;
exportToPDF = false;
renderVideo = false;

%% Create the optics
theOI = colorDetectOpticalImageConstruct(rParams.oiParams);

%% Create the cone mosaic
theMosaic = colorDetectConeMosaicConstruct(rParams.mosaicParams);

%% Define stimulus set
%
% Chromatic directions in L/M plane.  It's a little easier to think in
% terms of angles.  
%
% Then find highest within gamut vector for each direction
LMangles = (0:LMPlaneInstanceParams.deltaAngle:180-LMPlaneInstanceParams.deltaAngle)/180*pi;
for angleIndex = 1:numel(LMangles)
    theta = LMangles(angleIndex);
    baseTestConeContrastDirs(:,angleIndex) = LMPlaneInstanceParams.baseStimulusLength*[cos(theta) sin(theta) 0.0]';
    
    % Find the highest in gamut cone contrast and define cone contrast
    % vector to be just under this length.
    colorModulationParamsTemp = rParams.colorModulationParams;
    colorModulationParamsTemp.coneContrasts = baseTestConeContrastDirs(:,angleIndex);
    colorModulationParamsTemp.contrast = 1;
    [~,contrastScaleFactor(angleIndex)] = colorGaborSceneCreate(rParams.gaborParams,colorModulationParamsTemp,true);
    testConeContrasts(:,angleIndex) = 0.98*contrastScaleFactor(angleIndex)*baseTestConeContrastDirs(:,angleIndex);
end

% Contrasts
if (strcmp(LMPlaneInstanceParams.contrastScale, 'linear'))
    testContrasts = linspace(LMPlaneInstanceParams.lowContrast, LMPlaneInstanceParams.highContrast, LMPlaneInstanceParams.nContrastsPerDirection);
else
    testContrasts = logspace(log10(LMPlaneInstanceParams.lowContrast), log10(LMPlaneInstanceParams.highContrast), LMPlaneInstanceParams.nContrastsPerDirection);
end

%% Generate data for the no stimulus condition
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

% Write the no cone contrast data and some extra facts we need
currentParamsList = {rParams, colorModulationParamsTemp, LMPlaneInstanceParams};
rwObject.write('responseInstances',theNoStimData,parentParamsList,currentParamsList,theProgram);

%% Save the other data we need for use by the classifier preprocessing subroutine
ancillaryData = struct(...
    'testConeContrasts', testConeContrasts, ...
    'testContrasts', testContrasts, ...
    'theMosaic', theMosaic, ...
    'gaborParams', gaborParams, ...
    'temporalParams', temporalParams, ...
    'oiParams', oiParams, ...
    'mosaicParams', mosaicParams);
rwObject.write('ancillaryData',ancillaryData,parentParamsList,currentParamsList,theProgram);

%% Visualize responses

%% Generate data for all the examined stimuli
%
% This code is written in a slightly convoluted manner so that it will run
% inside a parfor loop -- we had to define some single indexed temp
% variables to get this to work.
%
% It is also possible that the parfor loop will not work for you, depending
% on your Matlab configuration.  In this case, change it to a for loop.

% Create crossed list of contrast directions and contrasts to run.
nParforConditions = size(testConeContrasts,2)*numel(testContrasts);
theParforConditionStructs = cell(nParforConditions,1);
conditionIndex = 1;
for ii = 1:size(testConeContrasts,2) 
    for jj = 1:numel(testContrasts)
        thisConditionStruct.ii = ii;
        thisConditionStruct.jj = jj;
        thisConditionStruct.testConeContrasts = testConeContrasts(:,ii);
        thisConditionStruct.contrast = testContrasts(:,ii);
        theParforConditionStructs{conditionIndex} = thisConditionStruct;
        conditionIndex = conditionIndex + 1;
    end
end

% Loop over color directions
tic;
parfor kk = 1:nParforConditions
    thisConditionStruct = theParforConditionStructs{kk};
    ii = thisConditionStruct.ii;
    jj = thisConditionStruct.jj;
    
    colorModulationParamsTemp = rParams.colorModulationParams;
    colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
    colorModulationParamsTemp.contrast = testContrasts;

    % Make noisy instances for each contrast
    stimDataJJ = cell(1,numel(testContrasts));
        stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.2f',...
            colorModulationParamsTemp.coneContrasts(1), colorModulationParamsTemp.coneContrasts(2), colorModulationParamsTemp.coneContrasts(3), colorModulationParamsTemp.contrast);
        stimData = struct(...
            'testContrast', colorModulationParamsTemp.contrast, ...
            'testConeContrasts', colorModulationParamsTemp.coneContrasts, ...
            'stimulusLabel', stimulusLabel, ...
            'responseInstanceArray', colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, LMPlaneInstanceParams.trialsNum, rParams.temporalParams.simulationTimeStepSecs, ...
            rParams.gaborParams, colorModulationParamsTemp, rParams.temporalParams, theOI, theMosaic)); 
        
    % Save data for this color direction/contrast pair
    currentParamsList = {rParams, colorModulationParamsTemp, LMPlaneInstanceParams};
    rwObject.write('responseInstances',theNoStimData,parentParamsList,currentParamsList,theProgram);
    %parforResponseInstancesSave(stimData,fullfile(outputDir,sprintf('responseInstances_%d_%d',thisConditionStruct.ii,thisConditionStruct.jj)));
end
fprintf('Finished generating responses in %2.2f minutes\n', toc/60);


% THIS IS BROKEN AND NEES TO BE UPDATED TO NEW DATA FORMAT AND rwObject
% land.
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
