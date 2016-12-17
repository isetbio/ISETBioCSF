function [validationData, extraData] = t_coneCurrentEyeMovementsResponseInstances(varargin)
% [validationData, extraData] = t_coneCurrentEyeMovementsResponseInstances(varargin)
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
%   'emPathType' - Value (one of: 'none', 'frozen', 'random') determines whether we have:
%                  zero eye movements across all trials, 
%                  an emPath that is frozen across all trials,
%                  or a dynamic emPath that changes across trials.
%   'freezeNoise' - true/false (default true).  Freezes all noise so that results are reproducible
%   'compute' - true/false (default true).  Do the computations.
%   'computeMosaic' - true/false (default true). Compute a cone mosaic or load one (good for large hex mosaics which take a while to compute)
%   'generatePlots' - true/false (default false).  Produce response
%        visualizations.  Set to false when running big jobs on clusters or
%        in parfor loops, as plotting doesn't seem to play well with those
%        conditions.
%   'visualizedResponseNormalization' - How to normalize visualized responses
%        Available options: 'submosaicBasedZscore', 'LMSabsoluteResponseBased', 'LMabsoluteResponseBased', 'MabsoluteResponseBased'
%   'exportPDF' - true/false (default true).  If visualizing responses,
%        export the PDF files.
%   'delete' - true/false (default true).  Delete the response instance
%        files.  Useful for cleaning up big output when we are done with
%        it.  If this is true, output files are deleted at the end.

%% Parse input
p = inputParser;
p.addParameter('rParams',[],@isemptyorstruct);
p.addParameter('testDirectionParams',[],@isemptyorstruct);
p.addParameter('emPathType','none',@ischar);
p.addParameter('freezeNoise',true,@islogical);
p.addParameter('compute',true,@islogical);
p.addParameter('computeMosaic', true, @islogical);
p.addParameter('generatePlots',false,@islogical);
p.addParameter('exportPDF',true,@islogical);
p.addParameter('delete',false',@islogical);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.addParameter('workerID', [], @isnumeric);
p.parse(varargin{:});
rParams = p.Results.rParams;
testDirectionParams = p.Results.testDirectionParams;

%% Clear
if (nargin == 0)
    ieInit; close all;
end

validationData = [];
extraData = [];

%% Get the parameters we need
%
% t_colorGaborResponseGenerationParams returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
if (isempty(rParams))
    rParams = responseParamsGenerate;
    
    % Override some defult parameters
    %
    % Set duration equal to sampling interval to do just one frame.
    rParams.temporalParams.stimulusDurationInSeconds = 200/1000;
    rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.stimulusDurationInSeconds;
    rParams.temporalParams.secondsToInclude = rParams.temporalParams.stimulusDurationInSeconds;
    
    rParams.mosaicParams.integrationTimeInSeconds = rParams.temporalParams.stimulusDurationInSeconds;
    rParams.mosaicParams.isomerizationNoise = 'random';         % Type coneMosaic.validNoiseFlags to get valid values
    rParams.mosaicParams.osNoise = 'random';                    % Type outerSegment.validNoiseFlags to get valid values
    rParams.mosaicParams.osModel = 'Linear';
end


% Set the emPathType
if (~isempty(p.Results.emPathType))
    rParams.temporalParams.emPathType = p.Results.emPathType;
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

%% Parameters that define the LM instances we'll generate here
if (isempty(testDirectionParams))
    testDirectionParams = instanceParamsGenerate;
end

% The constant params list
constantParamsList = {rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams,  rParams.backgroundParams, testDirectionParams};

colorModulationParamsNull = rParams.colorModulationParams;
colorModulationParamsNull.coneContrasts = [0 0 0]';
colorModulationParamsNull.contrast = 0;
    
%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
theProgram = mfilename;

%% The computing happens here, if we are doing it
if (p.Results.compute)
    % Create the optics
    theOI = colorDetectOpticalImageConstruct(rParams.oiParams);
    
    if (p.Results.computeMosaic)
        % Create the cone mosaic
        theMosaic = colorDetectConeMosaicConstruct(rParams.mosaicParams);
        % Save cone mosaic
        coneParamsList = {rParams.mosaicParams};
        rwObject.write('coneMosaic', theMosaic, coneParamsList, theProgram, 'type', 'mat');
    else
         % Load a previously saved cone mosaic
         fprintf(2,'Loading a previously saved cone mosaic\n');
         coneParamsList = {rParams.mosaicParams};
         theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    end
    
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
            
        case 'LMSPlane'
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
    fprintf('Computing the null stimulus response ...\n');
    stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.2f', ...
        colorModulationParamsNull.coneContrasts(1), colorModulationParamsNull.coneContrasts(2), colorModulationParamsNull.coneContrasts(3), colorModulationParamsNull.contrast);
    [responseInstanceArray,noiseFreeIsomerizations, noiseFreePhotocurrents] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, testDirectionParams.trialsNum, ...
        rParams.spatialParams, rParams.backgroundParams, colorModulationParamsNull, rParams.temporalParams, theOI, theMosaic, 'seed', 1, 'workerID', p.Results.workerID);
    
    noStimData = struct(...
        'testContrast', colorModulationParamsNull.contrast, ...
        'testConeContrasts', colorModulationParamsNull.coneContrasts, ...
        'stimulusLabel', stimulusLabel, ...
        'responseInstanceArray',responseInstanceArray, ...
        'noiseFreeIsomerizations',noiseFreeIsomerizations, ...
        'noiseFreePhotocurrents', noiseFreePhotocurrents);
    
    % Write the no cone contrast data and some extra facts we need
    paramsList = constantParamsList;
    paramsList{numel(paramsList)+1} = colorModulationParamsNull;
    rwObject.write('responseInstances',noStimData,paramsList,theProgram);
      
    % Save the other data we need for use by the classifier preprocessing subroutine
    ancillaryData = struct(...
        'testConeContrasts', testConeContrasts, ...
        'testContrasts', testContrasts, ...
        'rParams', rParams, ...
        'instanceParams', testDirectionParams);
    ancillaryData.parforConditionStructs = parforConditionStructs;
    rwObject.write('ancillaryData',ancillaryData,paramsList,theProgram);
    
    %% Generate data for all the examined stimuli
    %
    % It is possible that the parfor loop will not work for you, depending
    % on your Matlab configuration.  In this case, change it to a for loop.
   
    % Loop over color directions
    %
    % Note tha as the mosaic handle object (and any other handle objects)
    % enter the parfor loop, a copy local to each worker is created.  This
    % is the behavior we want, so that the isomerizations calculations for
    % each loop iteration don't step on each other.  Also note that any
    % changes to the mosaic object inside the loop do not get propagted
    % back out -- the mosaic object itself is the same at the end of the
    % parfor as at the start.  This is also OK here, but might be confusing
    % under other circumstances.
    tic;
    stimDataForValidation = cell(nParforConditions,1);

    parfor kk = 1:nParforConditions
        fprintf('Computing responses for condition %d/%d ...\n', kk,nParforConditions);
        if (~isempty(p.Results.workerID))
            % Get the parallel pool worker ID
            t = getCurrentTask();
            workerID = t.ID;
        else
            workerID = [];
        end
        
        thisConditionStruct = parforConditionStructs{kk};
        colorModulationParamsTemp = rParams.colorModulationParams;
        colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
        colorModulationParamsTemp.contrast = thisConditionStruct.contrast;
        
        % Make noisy instances for each contrast
        stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.5f',...
            colorModulationParamsTemp.coneContrasts(1), colorModulationParamsTemp.coneContrasts(2), colorModulationParamsTemp.coneContrasts(3), colorModulationParamsTemp.contrast);
        [responseInstanceArray,noiseFreeIsomerizations, noiseFreePhotocurrents] = ...
            colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, testDirectionParams.trialsNum, ...
                rParams.spatialParams, rParams.backgroundParams, colorModulationParamsTemp, rParams.temporalParams, theOI, theMosaic, 'seed', parforRanSeeds(kk), 'workerID', workerID);
        stimData = struct(...
            'testContrast', colorModulationParamsTemp.contrast, ...
            'testConeContrasts', colorModulationParamsTemp.coneContrasts, ...
            'stimulusLabel', stimulusLabel, ...
            'responseInstanceArray',responseInstanceArray, ...
            'noiseFreeIsomerizations',noiseFreeIsomerizations, ...
            'noiseFreePhotocurrents', noiseFreePhotocurrents);
        
        % Save some data for validation in first loop
        if (kk == 1)
            s = struct();
            savedTrial = 1;
            s.theMosaicIsomerizations = squeeze(responseInstanceArray.theMosaicIsomerizations(savedTrial,:,:,:));
            if (isempty(responseInstanceArray.theMosaicPhotocurrents))
                s.theMosaicPhotocurrents = [];
            else
                s.theMosaicPhotocurrents = squeeze(responseInstanceArray.theMosaicPhotocurrents(savedTrial,:,:,:));
            end
            s.theMosaicEyeMovements = squeeze(responseInstanceArray.theMosaicEyeMovements(savedTrial,:,:));
            s.timeAxis = responseInstanceArray.timeAxis;
            s.photocurrentTimeAxis = responseInstanceArray.photocurrentTimeAxis;
            stimDataForValidation{kk} = s;
        end
        
        % Save data for this color direction/contrast pair
        paramsList = constantParamsList;
        paramsList{numel(paramsList)+1} = colorModulationParamsTemp;
        rwObject.write('responseInstances',stimData,paramsList,theProgram);
    end
    fprintf('Finished generating responses in %2.2f minutes\n', toc/60);
    
    %% Validation data
    %
    % Since we want to store some stimData and that isn't saved outside the
    % parfor loop (and is in fact a bit hard to save outside the parfor
    % loop)
    if (nargout > 0)
        savedTrial = 1;
        noStimValidationData.theMosaicIsomerizations = squeeze(noStimData.responseInstanceArray.theMosaicIsomerizations(savedTrial,:,:,:));
        if (isempty(noStimData.responseInstanceArray.theMosaicPhotocurrents))
            noStimValidationData.theMosaicPhotocurrents = [];
        else
            noStimValidationData.theMosaicPhotocurrents = squeeze(noStimData.responseInstanceArray.theMosaicPhotocurrents(savedTrial,:,:,:));
        end
        noStimValidationData.theMosaicEyeMovements = squeeze(noStimData.responseInstanceArray.theMosaicEyeMovements(savedTrial,:,:));
        noStimValidationData.timeAxis = noStimData.responseInstanceArray.timeAxis;
        noStimValidationData.photocurrentTimeAxis = noStimData.responseInstanceArray.photocurrentTimeAxis;
        
        validationData.noStimData = noStimValidationData;
        validationData.stimData = stimDataForValidation;
        extraData.ancillaryData = ancillaryData;
        extraData.p.Results = p.Results;
    end
end

%% Visualize
if (p.Results.generatePlots)

    % How many istances to visualize
    instancesToVisualize = 1:5;
    
    % Load the mosaic
    coneParamsList = {rParams.mosaicParams};
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
         
    % Load the response and ancillary data
    paramsList = constantParamsList;
    paramsList{numel(paramsList)+1} = colorModulationParamsNull;    
    noStimData = rwObject.read('responseInstances',paramsList,theProgram);
    ancillaryData = rwObject.read('ancillaryData',paramsList,theProgram);
    
    rParams = ancillaryData.rParams;
    parforConditionStructs = ancillaryData.parforConditionStructs;
    nParforConditions = length(parforConditionStructs); 
    for kk = 1:nParforConditions 
         thisConditionStruct = parforConditionStructs{kk};
         colorModulationParamsTemp = rParams.colorModulationParams;
         colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
         colorModulationParamsTemp.contrast = thisConditionStruct.contrast;

         paramsList = constantParamsList;
         paramsList{numel(paramsList)+1} = colorModulationParamsTemp;    
         stimData = rwObject.read('responseInstances',paramsList,theProgram);
         visualizeResponseInstances(theMosaic, stimData, noStimData, p.Results.visualizedResponseNormalization, kk, nParforConditions, instancesToVisualize);
    end
end


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
        paramsList = constantParamsList;
        paramsList{numel(paramsList)+1} = colorModulationParamsTemp;
    
       % paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
        rwObject.delete('responseInstances',paramsList,theProgram);   
    end   
end
end
