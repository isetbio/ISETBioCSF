function validationData = t_colorGaborConeCurrentEyeMovementsResponseInstances(rParams,testDirectionParams)
% validationData = t_colorGaborConeCurrentEyeMovementsResponseInstances([rParams],[testDirectionParams])
%
% Show how to generate a number of response instances for a given stimulus
% condition.  The default parameters are set up to generate just a single frame
% of the response, but the same tutorial can do temporal sequences with other 
% parameter choices.
%
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
if (nargin == 0)
    ieInit; close all;
end

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
    % Set duration equal to sampling interval to do just one frame.
    rParams.temporalParams.simulationTimeStepSecs = 200/1000;
    rParams.temporalParams.stimulusDurationInSeconds = rParams.temporalParams.simulationTimeStepSecs;
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
if (nargin < 2 | isempty(testDirectionParams))
    testDirectionParams = LMPlaneInstanceParamsGenerate;
end

%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
theProgram = mfilename;

%% Control what types of output are written
% These may only work on some computers, depending on what infrastructure is installed.
visualizeResponses = false;
exportToPDF = false;
renderVideo = false;

%% Create the optics
theOI = colorDetectOpticalImageConstruct(rParams.oiParams);

%% Create the cone mosaic
rParams.mosaicParams.fieldOfViewDegs = rParams.gaborParams.fieldOfViewDegs;
theMosaic = colorDetectConeMosaicConstruct(rParams.mosaicParams);

%% Define color direction cone contrasts as well as contrast scalars.
%
% Directions
testConeContrasts = testConeContrastsFromTestDirectionParams(rParams,testDirectionParams);

% Contrasts
if (strcmp(testDirectionParams.contrastScale, 'linear'))
    testContrasts = linspace(testDirectionParams.lowContrast, testDirectionParams.highContrast, testDirectionParams.nContrastsPerDirection);
else
    testContrasts = logspace(log10(testDirectionParams.lowContrast), log10(testDirectionParams.highContrast), testDirectionParams.nContrastsPerDirection);
end

%% Generate data for the no stimulus condition
colorModulationParamsTemp = rParams.colorModulationParams;
colorModulationParamsTemp.coneContrasts = [0 0 0]';
colorModulationParamsTemp.contrast = 0;
stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.2f', ...
    colorModulationParamsTemp.coneContrasts(1), colorModulationParamsTemp.coneContrasts(2), colorModulationParamsTemp.coneContrasts(3), colorModulationParamsTemp.contrast);
[responseInstanceArray,noiseFreeIsomerizations] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, testDirectionParams.trialsNum, rParams.temporalParams.simulationTimeStepSecs, ...
    rParams.gaborParams, rParams.backgroundParams, colorModulationParamsTemp, rParams.temporalParams, theOI, theMosaic);
noStimData = struct(...
    'testContrast', colorModulationParamsTemp.contrast, ...
    'testConeContrasts', colorModulationParamsTemp.coneContrasts, ...
    'stimulusLabel', stimulusLabel, ...
    'responseInstanceArray',responseInstanceArray, ...
    'noiseFreeIsomerizations',noiseFreeIsomerizations);

% Write the no cone contrast data and some extra facts we need
paramsList = {rParams.gaborParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
rwObject.write('responseInstances',noStimData,paramsList,theProgram);

%% Save the other data we need for use by the classifier preprocessing subroutine
ancillaryData = struct(...
    'testConeContrasts', testConeContrasts, ...
    'testContrasts', testContrasts, ...
    'theMosaic', theMosaic, ...
    'rParams', rParams, ...
    'LMPlaneInstanceParams', testDirectionParams);
rwObject.write('ancillaryData',ancillaryData,paramsList,theProgram);

%% Generate data for all the examined stimuli
%
% It is possible that the parfor loop will not work for you, depending
% on your Matlab configuration.  In this case, change it to a for loop.
 
% Loop over color directions
tic;
parforConditionStructs = responseGenerationParforConditionStructsGenerate(testConeContrasts,testContrasts);
nParforConditions = length(parforConditionStructs);
parfor kk = 1:nParforConditions
    thisConditionStruct = parforConditionStructs{kk};
    colorModulationParamsTemp = rParams.colorModulationParams;
    colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
    colorModulationParamsTemp.contrast = thisConditionStruct.contrast;
    
    % Make noisy instances for each contrast
    stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.2f',...
        colorModulationParamsTemp.coneContrasts(1), colorModulationParamsTemp.coneContrasts(2), colorModulationParamsTemp.coneContrasts(3), colorModulationParamsTemp.contrast);
    [responseInstanceArray,noiseFreeIsomerizations] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, testDirectionParams.trialsNum, rParams.temporalParams.simulationTimeStepSecs, ...
        rParams.gaborParams, rParams.backgroundParams, colorModulationParamsTemp, rParams.temporalParams, theOI, theMosaic);
    stimData = struct(...
        'testContrast', colorModulationParamsTemp.contrast, ...
        'testConeContrasts', colorModulationParamsTemp.coneContrasts, ...
        'stimulusLabel', stimulusLabel, ...
        'responseInstanceArray',responseInstanceArray, ...
        'noiseFreeIsomerizations',noiseFreeIsomerizations);
    
    % Save data for this color direction/contrast pair
    paramsList = {rParams.gaborParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
    rwObject.write('responseInstances',stimData,paramsList,theProgram);
end
fprintf('Finished generating responses in %2.2f minutes\n', toc/60);

%% Validation data
if (nargin > 0)
    validationData = [];
end

% THIS IS BROKEN AND NEES TO BE UPDATED TO NEW DATA FORMAT AND rwObject
% LAND.
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
                figHandle = visualizeResponseInstance(conditionDir, s.responseInstanceArray(iTrial), stimulusLabel, theMosaic, iTrial, testDirectionParams.trialsNum, renderVideo);
                if (exportToPDF)
                    figFileNames{ii, jj, iTrial} = ...
                        fullfile(colorGaborDetectOutputDir(conditionDir),sprintf('%s_Trial%dOf%d.pdf', stimulusLabel, iTrial, testDirectionParams.trialsNum),'figures');
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
