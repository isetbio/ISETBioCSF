function validationData = t_coneCurrentEyeMovementsResponseInstancesOLD(varargin)
% validationData = t_coneCurrentEyeMovementsResponseInstances(varargin)
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
%   t_coneCurrentEyeMovementsMovie
% but the actual routine has some tricks to make it go fast.  There is also
% are routine
%   colorDetectResponseInstanceArrayConstruct
% that works more like the tutorial but is slower.
%
% This tutorial saves its output in a .mat file, which cah then read in by
%   t_colorDetectFindPerformance
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
%   t_coneCurrentEyeMovementsMovie
%	t_colorDetectFindPerformance
%   colorDetectResponseInstanceArrayConstruct
%   colorDetectResponseInstanceFastArrayConstruct
%
% Key/value pairs
%   'rParams' - Value the is the rParams structure to use
%   'testDirectionParams - Value is the testDirectionParams structure to use
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
p.addParameter('testDirectionParams',[],@isemptyorstruct);
p.addParameter('setRng',true,@islogical);
p.addParameter('compute',true,@islogical);
p.addParameter('generatePlots',false,@islogical);
p.addParameter('exportPDF',true,@islogical);
p.addParameter('renderVideo',true,@islogical);
p.addParameter('delete',false',@islogical);
p.parse(varargin{:});
rParams = p.Results.rParams;
testDirectionParams = p.Results.testDirectionParams;

%% Clear
if (nargin == 0)
    ieInit; close all;
end

%% Fix random number generator so we can validate output exactly
if (p.Results.setRng)
    rng(1);
end

%% Get the parameters we need
%
% t_colorGaborResponseGenerationParams returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
if (isempty(rParams))
    rParams = responseParamsGenerate;
    
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
if (isempty(testDirectionParams))
    testDirectionParams = instanceParamsGenerate;
end

%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
theProgram = mfilename;

%% The computing happens here, if we are doing it
if (p.Results.compute)
    % Create the optics
    theOI = colorDetectOpticalImageConstruct(rParams.oiParams);
    
    % Create the cone mosaic
    rParams.mosaicParams.fieldOfViewDegs = rParams.spatialParams.fieldOfViewDegs;
    theMosaic = colorDetectConeMosaicConstruct(rParams.mosaicParams);
    
    %% Define color modulation list
    switch (testDirectionParams.instanceType)
        case 'LMPlane'
            % Directions
            testConeContrasts = testConeContrastsFromTestDirectionParams(rParams,testDirectionParams);
            
            % Contrasts
            if (strcmp(testDirectionParams.contrastScale, 'linear'))
                testContrasts = linspace(testDirectionParams.lowContrast, testDirectionParams.highContrast, testDirectionParams.nContrastsPerDirection);
            else
                testContrasts = logspace(log10(testDirectionParams.lowContrast), log10(testDirectionParams.highContrast), testDirectionParams.nContrastsPerDirection);
            end
            
        case 'contrasts'
            % Contrasts
            if (strcmp(testDirectionParams.contrastScale, 'linear'))
                testContrasts = linspace(testDirectionParams.lowContrast, testDirectionParams.highContrast, testDirectionParams.nContrastsPerDirection);
            else
                testContrasts = logspace(log10(testDirectionParams.lowContrast), log10(testDirectionParams.highContrast), testDirectionParams.nContrastsPerDirection);
            end
            
            % Set up parfor condition structs. Just dummy up a cone
            % contrasts argument for this call, we won't use the returned
            % contrasts for this case below.
            testConeContrasts = NaN*zeros(3,1);
        otherwise
            error('Unknown test instance type passed');
    end
    parforConditionStructs = responseGenerationParforConditionStructsGenerate(testConeContrasts,testContrasts);
    nParforConditions = length(parforConditionStructs);
    parforRanSeeds = randi(1000000,nParforConditions,1)+1;
    
    % Generate data for the no stimulus condition
    colorModulationParamsTemp = rParams.colorModulationParams;
    colorModulationParamsTemp.coneContrasts = [0 0 0]';
    colorModulationParamsTemp.contrast = 0;
    stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.2f', ...
        colorModulationParamsTemp.coneContrasts(1), colorModulationParamsTemp.coneContrasts(2), colorModulationParamsTemp.coneContrasts(3), colorModulationParamsTemp.contrast);
    [responseInstanceArray,noiseFreeIsomerizations] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, testDirectionParams.trialsNum, rParams.temporalParams.simulationTimeStepSecs, ...
        rParams.spatialParams, rParams.backgroundParams, colorModulationParamsTemp, rParams.temporalParams, theOI, theMosaic);
    noStimData = struct(...
        'testContrast', colorModulationParamsTemp.contrast, ...
        'testConeContrasts', colorModulationParamsTemp.coneContrasts, ...
        'stimulusLabel', stimulusLabel, ...
        'responseInstanceArray',responseInstanceArray, ...
        'noiseFreeIsomerizations',noiseFreeIsomerizations);
    
    % Write the no cone contrast data and some extra facts we need
    paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
    rwObject.write('responseInstances',noStimData,paramsList,theProgram);
    
    % Save the other data we need for use by the classifier preprocessing subroutine
    ancillaryData = struct(...
        'testConeContrasts', testConeContrasts, ...
        'testContrasts', testContrasts, ...
        'rParams', rParams, ...
        'instanceParams', testDirectionParams);
    rwObject.write('ancillaryData',ancillaryData,paramsList,theProgram);
    
    %% Generate data for all the examined stimuli
    %
    % It is possible that the parfor loop will not work for you, depending
    % on your Matlab configuration.  In this case, change it to a for loop.

    % Loop over color directions
    tic;
    stimDataForValidation = cell(nParforConditions,1);
    rState = rng;
    parfor kk = 1:nParforConditions
        rng(parforRanSeeds(kk));
        thisConditionStruct = parforConditionStructs{kk};
        colorModulationParamsTemp = rParams.colorModulationParams;
        colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
        colorModulationParamsTemp.contrast = thisConditionStruct.contrast;
        
        % Make noisy instances for each contrast
        stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.2f',...
            colorModulationParamsTemp.coneContrasts(1), colorModulationParamsTemp.coneContrasts(2), colorModulationParamsTemp.coneContrasts(3), colorModulationParamsTemp.contrast);
        [responseInstanceArray,noiseFreeIsomerizations] = colorDetectResponseInstanceArrayFastConstructOLD(stimulusLabel, testDirectionParams.trialsNum, rParams.temporalParams.simulationTimeStepSecs, ...
            rParams.spatialParams, rParams.backgroundParams, colorModulationParamsTemp, rParams.temporalParams, theOI, theMosaic);
        stimData = struct(...
            'testContrast', colorModulationParamsTemp.contrast, ...
            'testConeContrasts', colorModulationParamsTemp.coneContrasts, ...
            'stimulusLabel', stimulusLabel, ...
            'responseInstanceArray',responseInstanceArray, ...
            'noiseFreeIsomerizations',noiseFreeIsomerizations);
        
        % Save some data for validation in first loop
        if (kk == 1)
            stimDataForValidation{kk} = stimData.responseInstanceArray(1);
        end
        
        % Save data for this color direction/contrast pair
        paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
        rwObject.write('responseInstances',stimData,paramsList,theProgram);
    end
    rng(rState);
    fprintf('Finished generating responses in %2.2f minutes\n', toc/60);
    
    %% Validation data
    %
    % Since we want to store some stimData and that isn't saved outside the
    % parfor loop (and is in fact a bit hard to save outside the parfor
    % loop)
    if (nargout > 0)
        validationData.ancillaryData = ancillaryData;
        validationData.noStimData = noStimData.responseInstanceArray(1);
        validationData.stimData = stimDataForValidation{1};
    end
end

%% Visualize
%
% This is not yet implemented, but here is where the code would go.  And, a
% very early version of the visualization is in the commented out code
% immediately below
if (p.Results.generatePlots)
    % noStimData = rwObject.read('responseInstances',paramsList,theProgram);
    % ancillaryData = rwObject.read('ancillaryData',paramsList,theProgram);
    % parforConditionStructs = ancillaryData.parforConditionStructs;
    % nParforConditions = length(parforConditionStructs);
    % 
    % for kk = 1:nParforConditions 
    %     thisConditionStruct = parforConditionStructs{kk};
    %     colorModulationParamsTemp = rParams.colorModulationParams;
    %     colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
    %     colorModulationParamsTemp.contrast = thisConditionStruct.contrast;
    %     paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
    %     stimData = rwObject.read('responseInstances',paramsList,theProgram);   
    % end   
end
% % THIS IS BROKEN AND NEES TO BE UPDATED TO NEW DATA FORMAT AND rwObject
% % LAND.
% %
% % Also, the time numbers on the videos do not seem to correspond to the
% % stimulus peak at time 0, and the videos look screwy in the no noise
% % condition.
% if (p.Results.generatePlots)
    % fprintf('\nVisualizing responses ...\n');
    % for ii = 1:size(testConeContrasts,2)
    %     for jj = 1:numel(testContrasts)
    %         stimulusLabel = sprintf('LMS_%2.2f_%2.2f_%2.2f_Contrast_%2.2f', testConeContrasts(1,ii), testConeContrasts(2,ii), testConeContrasts(3,ii), testContrasts(jj));
    %         s = theStimData{ii, jj};
    % 
    %         % Visualize a few response instances only
    %         for iTrial = 1:2
    %             figHandle = visualizeResponseInstance(conditionDir, s.responseInstanceArray(iTrial), stimulusLabel, theMosaic, iTrial, testDirectionParams.trialsNum, p.Results.renderVideo);
    %             if (p.Results.exportPDF)
    %                 figFileNames{ii, jj, iTrial} = ...
    %                     fullfile(colorGaborDetectOutputDir(conditionDir),sprintf('%s_Trial%dOf%d.pdf', stimulusLabel, iTrial, testDirectionParams.trialsNum),'figures');
    %                 NicePlot.exportFigToPDF(figFileNames{ii, jj, iTrial}, figHandle, 300);
    %             end
    %         end
    %     end
    % end
    % 
    % % Export summary PDF with all responses
    % if (p.Results.exportPDF)
    %     summaryPDF = fullfile(colorGaborDetectOutputDir(conditionDir,'figures'), 'AllInstances.pdf');
    %     fprintf('Exporting a summary PDF with all response instances in %s\n', summaryPDF);
    %     NicePlot.combinePDFfilesInSinglePDF(figFileNames(:), summaryPDF);
    % end
% end

%% Delete response instance files.
%
% Doesn't delete figures or parent directories.
if (p.Results.delete)
    rwObject.delete('responseInstances',paramsList,theProgram);
    rwObject.delete('ancillaryData',paramsList,theProgram);
    for kk = 1:nParforConditions 
        thisConditionStruct = parforConditionStructs{kk};
        colorModulationParamsTemp = rParams.colorModulationParams;
        colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
        colorModulationParamsTemp.contrast = thisConditionStruct.contrast;
        paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
        rwObject.delete('responseInstances',paramsList,theProgram);   
    end   
end
