function [validationData, extraData] = t_coneCurrentEyeMovementsResponseInstances(varargin)
% T_CONECURRENTEYEMOVEMENTRESPONSES  Generate response instances for a given stimulus condition.
%     [validationData, extraData] = T_CONECURRENTEYEMOVEMENTRESPONSES(varargin)
% 
%     Show how to generate a number of response instances for a given
%     stimulus condition.  The default parameters are set up to generate
%     just a single frame of the response, but the same tutorial can do
%     temporal sequences with other parameter choices.
% 
%     This tutorial relies on routine
%     colorDetectResponseInstanceFastArrayConstruct which does most of the
%     hard work.  The basic principles underlying
%     colorDetectResponseInstanceFastArrayConstruct itself is demonstrated
%     in tutorial t_coneCurrentEyeMovementsMovie but the actual routine has
%     some tricks to make it go fast.  There is also are routine
%     colorDetectResponseInstanceArrayConstruct that works more like the
%     tutorial but is slower.
% 
%     This tutorial saves its output in a .mat file, which cah then read in
%     by t_colorDetectFindPerformance which shows how to use the data to
%     find the thresholds.
% 
%     The returned validation structure allows this routine to be called
%     from a validation script driven by the UnitTest toolbox.
% 
%     The tutorial produces output according to a scheme controlled by the
%     specified IBIOColorDetect rwObject.
% 
%     Key/value pairs
%     'rParams' - Value the is the rParams structure to use.  Default empty,
%     which then uses defaults produced by generation function.
%     'testDirectionParams - Value is the testDirectionParams structure to use
%     'centeredEMPaths' - true/false (default false) 
%               Controls wether the eye movement paths start at (0,0) (default) or wether they are centered around (0,0)
%     'trialBlockSize' - How many blocks to split the testDirectionParams.trialsNum into. Default: [], which results in nTrials  (no blocking). 
%               This only has an effect with @coneMosaicHex mosaics and when nTrials>1 and it is useful with 
%               large mosaics x lots of trials, in which case the absorptions matrix does not fit in the RAM.
%               If set to -1, the number of trial blocks is computed automatically based on the number of cores and system RAM.
%     'displayTrialBlockPartitionDiagnostics', true/false. Wether to display trial block diagnostics.
%     'freezeNoise' - true/false (default true).  Freezes all noise so that results are reproducible
%     'compute' - true/false (default true).  Do the computations.
%     'computeMosaic' - true/false (default true). Compute a cone mosaic or load one (good for large hex mosaics which take a while to compute)
%     'parforWorkersNum' - 0 .. 12 (default: 12). How many workers to use for the computations.
%       use 0: for a serial for loop
%       use > 0: for a parfor loop with desired number of workers
%     'generatePlots' - true/false (default false).  Produce response
%        visualizations.  Set to false when running big jobs on clusters or
%        in parfor loops, as plotting doesn't seem to play well with those
%        conditions.
%     'visualizeResponses' - true/false (default true). Call the fancy visualize response routine?
%     'visualizationFormat' - How to arrange visualized maps. 
%       Available options: 'montage', 'video'. Default is 'montage'
%     'visualizedResponseNormalization' - How to normalize visualized responses
%        Available options: 'submosaicBasedZscore', 'LMSabsoluteResponseBased', 'LMabsoluteResponseBased', 'MabsoluteResponseBased'
%     'exportPDF' - true/false (default true).  If visualizing responses,
%        export the PDF files.
%     'delete' - true/false (default true).  Delete the response instance
%        files.  Useful for cleaning up big output when we are done with
%        it.  If this is true, output files are deleted at the end.
%
%     See also: T_COLORDETECTFINDPERFORMANCE COLORDETECTRESPONSEINSTANCEFASTARRAYCONSTRUCT

%% Parse input
p = inputParser;
p.addParameter('rParams',[],@isemptyorstruct);
p.addParameter('testDirectionParams',[],@isemptyorstruct);
p.addParameter('centeredEMPaths',false, @islogical); 
p.addParameter('trialBlockSize', [], @isnumeric);
p.addParameter('displayTrialBlockPartitionDiagnostics', false, @islogical);
p.addParameter('freezeNoise',true,@islogical);
p.addParameter('compute',true,@islogical);
p.addParameter('computeMosaic', true, @islogical);
p.addParameter('parforWorkersNum', 12, @isnumeric);
p.addParameter('overrideMosaicIntegrationTime', [], @isnumeric);
p.addParameter('generatePlots',false,@islogical);
p.addParameter('visualizeResponses',true,@islogical);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.addParameter('visualizationFormat', 'montage', @ischar);
p.addParameter('workerID', [], @isnumeric);
p.addParameter('exportPDF',true,@islogical);
p.addParameter('delete',false',@islogical);
p.parse(varargin{:});
rParams = p.Results.rParams;
testDirectionParams = p.Results.testDirectionParams;
visualizationFormat = p.Results.visualizationFormat;

% Ensure visualizationFormat has a valid value
if (strcmp(visualizationFormat, 'montage')) || (strcmp(visualizationFormat, 'video'))
else
    error('visualizationFormat must be set to either ''montage'' or ''video''. Current value: ''%s''.', visualizationFormat);
end

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
    
%% Fix random number generator so we can validate output exactly
if (p.Results.freezeNoise)
     fprintf(1, '\n%s: freezing all noise \n', mfilename);
     rng(1);
     if (strcmp(rParams.mosaicParams.isomerizationNoise, 'random'))
         fprintf(1, '\tmosaicParams.isomerizationNoise was set to ''%s'', setting it to ''frozen''.\n', rParams.mosaicParams.isomerizationNoise);
         rParams.mosaicParams.isomerizationNoise = 'frozen';
     end
     if (strcmp(rParams.mosaicParams.osNoise, 'random'))
         fprintf(1, '\tmosaicParams.osNoise was set to ''%s'', setting it to ''frozen''.\n', rParams.mosaicParams.osNoise);
         rParams.mosaicParams.osNoise = 'frozen';
     end
end

%% Parameters that define the LM instances we'll generate here
if (isempty(testDirectionParams))
    testDirectionParams = instanceParamsGenerate;
end

%% The constant params list
constantParamsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams,  rParams.backgroundParams, testDirectionParams};

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
        coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
        rwObject.write('coneMosaic', theMosaic, coneParamsList, theProgram, 'type', 'mat');
    else
         % Load a previously saved cone mosaic
         fprintf('Loading a previously saved cone mosaic\n');
         coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
         theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    end
    
    if (~isempty(p.Results.overrideMosaicIntegrationTime))
        fprintf(2, 'NOTE: Overriding mosaic''s default integrationTime (%2.1fms) with %2.1fms\n', theMosaic.integrationTime*1000, p.Results.overrideMosaicIntegrationTime*1000);
        theMosaic.integrationTime = p.Results.overrideMosaicIntegrationTime;
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
        rParams.spatialParams, rParams.backgroundParams, colorModulationParamsNull, rParams.temporalParams, theOI, theMosaic, ...
        'centeredEMPaths', p.Results.centeredEMPaths, ...
        'seed', 1, ...
        'workerID', p.Results.workerID,...
        'displayTrialBlockPartitionDiagnostics', p.Results.displayTrialBlockPartitionDiagnostics, ...
        'trialBlockSize', []);   % use all trials since this is done outside the parfor loop
    
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
        
    parfor (kk = 1:nParforConditions, p.Results.parforWorkersNum)
        fprintf('Computing responses for condition %d/%d ...\n', kk,nParforConditions);
        
        % Get the parallel pool worker ID
        if (~isempty(p.Results.workerID))  
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
                rParams.spatialParams, rParams.backgroundParams, colorModulationParamsTemp, rParams.temporalParams, theOI, theMosaic, ...
                'centeredEMPaths', p.Results.centeredEMPaths, ...
                'seed', parforRanSeeds(kk), ...
                'workerID', workerID, ...
                'displayTrialBlockPartitionDiagnostics', false, ...
                'trialBlockSize', p.Results.trialBlockSize);
            
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
if (p.Results.generatePlots && p.Results.visualizeResponses)

    % How many istances to visualize
    instancesToVisualize = 5;
    
    % Load the mosaic
    coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
         
    % Load the response and ancillary data
    paramsList = constantParamsList;
    paramsList{numel(paramsList)+1} = colorModulationParamsNull;   
    fprintf('Importing NOSTIM data and ancillary data\n');
    noStimData = rwObject.read('responseInstances',paramsList,theProgram);
    ancillaryData = rwObject.read('ancillaryData',paramsList,theProgram);
    
    % Only keep the data we will visualize
    noStimData.responseInstanceArray.theMosaicIsomerizations = noStimData.responseInstanceArray.theMosaicIsomerizations(1:instancesToVisualize,:,:);
    if (~isempty(noStimData.responseInstanceArray.theMosaicPhotocurrents))
        noStimData.responseInstanceArray.theMosaicPhotocurrents  = noStimData.responseInstanceArray.theMosaicPhotocurrents(1:instancesToVisualize,:,:);
    end
    noStimData.responseInstanceArray.theMosaicEyeMovements   = noStimData.responseInstanceArray.theMosaicEyeMovements(1:instancesToVisualize,:,:);
    rParams = ancillaryData.rParams;
    parforConditionStructs = ancillaryData.parforConditionStructs;
    nParforConditions = length(parforConditionStructs); 
    for kk = nParforConditions:-nParforConditions+1:1
        fprintf('Importing STIM data for condition %d/%d\n', kk, nParforConditions);
        thisConditionStruct = parforConditionStructs{kk};
        colorModulationParamsTemp = rParams.colorModulationParams;
        colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
        colorModulationParamsTemp.contrast = thisConditionStruct.contrast;

        paramsList = constantParamsList;
        paramsList{numel(paramsList)+1} = colorModulationParamsTemp;    
        stimData = rwObject.read('responseInstances',paramsList,theProgram);
        % Only keep the data we will visualize
        stimData.responseInstanceArray.theMosaicIsomerizations = stimData.responseInstanceArray.theMosaicIsomerizations(1:instancesToVisualize,:,:);
        if (~isempty(stimData.responseInstanceArray.theMosaicPhotocurrents))
            stimData.responseInstanceArray.theMosaicPhotocurrents  = stimData.responseInstanceArray.theMosaicPhotocurrents(1:instancesToVisualize,:,:);
        end
        stimData.responseInstanceArray.theMosaicEyeMovements   = stimData.responseInstanceArray.theMosaicEyeMovements(1:instancesToVisualize,:,:);
        visualizeResponseInstances(theMosaic, stimData, noStimData, p.Results.visualizedResponseNormalization, kk, nParforConditions, p.Results.visualizationFormat);
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
