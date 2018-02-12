function [validationData, extraData, varargout] = t_coneCurrentEyeMovementsResponseInstances(varargin)
% T_CONECURRENTEYEMOVEMENTRESPONSES  Generate response instances for a given stimulus condition.
%     [validationData, extraData, varargout] = T_CONECURRENTEYEMOVEMENTRESPONSES(varargin)
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
%     'displayTrialBlockPartitionDiagnostics', true/false. Wether to display trial block diagnostics.
%     'freezeNoise' - true/false (default true).  Freezes all noise so that results are reproducible
%     'compute' - true/false (default true).  Do the computations.
%     'computeOptics' - true/false (default true). Compute the optics
%     'computePhotocurrentResponseInstances' - true/false (default true) wether to compute photocurrent response instances
%     'computeMosaic' - true/false (default true). Compute a cone mosaic or load one (good for large hex mosaics which take a while to compute)
%     'ramPercentageEmployed' [>0 .. <=1]. Percentage of RAM to use for the computations. 
%       Use 0.5, if you want to run 2 simultaneous instances of MATLAB (touse all cores for example).
%       Use 1.0, otherwise
%     'parforWorkersNum' - 0 .. 20 (default: 20). How many workers to use for the computations.
%       use 0: for a serial for loop
%       use > 0: for a parfor loop with desired number of workers
%     'employStandardHostComputerResources' - true/false (default false). Only validation scripts should set this to true, so as to produce
%               identical sequences of random numbers. All other scripts should set it (leave it) to false.
%     'generatePlots' - true/false (default false).  Produce response
%        visualizations.  Set to false when running big jobs on clusters or
%        in parfor loops, as plotting doesn't seem to play well with those
%        conditions.
%     'visualizeMosaic' - true/false (default true). Wether to visualize the cone mosaic
%     'visualizeOptics' - true/false (default true). Wether to visualize the cone mosaic
%     'visualizeSpatialScheme' - true/false (default false). Visualize the relationship between mosaic and stimulus.
%     'visualizeOIsequence' - true/false (default false). Visualize the sequence of optical images
%     'visualizeResponses' - true/false (default true). Call the fancy visualize response routine
%     'visualizedConditionIndices' - list of conditions indices for which to visualize respondrd (default: empty - all conditions)
%     'visualizeOuterSegmentFilters' - true/false (default false). Visualize the outer segment impulse response functions.
%     'visualizationFormat' - How to arrange visualized maps. 
%       Available options: 'montage', 'video'. Default is 'montage'
%     'visualizedResponseNormalization' - How to normalize visualized responses
%        Available options: 'submosaicBasedZscore', 'LMSabsoluteResponseBased', 'LMabsoluteResponseBased', 'MabsoluteResponseBased'
%     'exportPDF' - true/false (default true).  If visualizing responses,
%        export the PDF files.
%     'delete' - true/false (default false).  Delete the response instance
%        files.  Useful for cleaning up big output when we are done with
%        it.  If this is true, output files are deleted at the end.
%
%     See also: T_COLORDETECTFINDPERFORMANCE COLORDETECTRESPONSEINSTANCEFASTARRAYCONSTRUCT

%% Parse input
p = inputParser;
p.addParameter('rParams',[],@isemptyorstruct);
p.addParameter('testDirectionParams',[],@isemptyorstruct);
p.addParameter('centeredEMPaths',false, @islogical); 
p.addParameter('displayTrialBlockPartitionDiagnostics', false, @islogical);
p.addParameter('freezeNoise',true,@islogical);
p.addParameter('compute',true,@islogical);
p.addParameter('computePhotocurrentResponseInstances', true, @islogical);
p.addParameter('computeMosaic', true, @islogical);
p.addParameter('visualizeMosaic',true, @islogical);
p.addParameter('computeOptics', true, @islogical);
p.addParameter('visualizeOptics',true, @islogical);
p.addParameter('ramPercentageEmployed', 1.0, @isnumeric);
p.addParameter('parforWorkersNum', 20, @isnumeric);
p.addParameter('employStandardHostComputerResources', false, @islogical);
p.addParameter('overrideMosaicIntegrationTime', [], @isnumeric);
p.addParameter('generatePlots',false,@islogical);
p.addParameter('visualizeResponses',true,@islogical);
p.addParameter('visualizedConditionIndices', [], @isnumeric);
p.addParameter('visualizeOuterSegmentFilters',false, @islogical);
p.addParameter('visualizeSpatialScheme',false,@islogical);
p.addParameter('visualizeOIsequence', false, @islogical);
p.addParameter('visualizeMosaicWithFirstEMpath', false, @islogical);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.addParameter('visualizationFormat', 'montage', @ischar);
p.addParameter('workerID', [], @isnumeric);
p.addParameter('exportPDF',true,@islogical);
p.addParameter('delete',false,@islogical);
p.addParameter('IBIOColorDetectSnapshot', [], @isstruct);

p.parse(varargin{:});
rParams = p.Results.rParams;
testDirectionParams = p.Results.testDirectionParams;
visualizationFormat = p.Results.visualizationFormat;
parforWorkersNum = p.Results.parforWorkersNum;

% Check the parforWorkersNum
[numberOfWorkers, ~, ~] = determineSystemResources(p.Results.employStandardHostComputerResources);
if (numberOfWorkers < parforWorkersNum)
    parforWorkersNum = numberOfWorkers;
end

% Ensure visualizationFormat has a valid value
if (strcmp(visualizationFormat, 'montage')) || (strcmp(visualizationFormat, 'video'))
else
    error('visualizationFormat must be set to either ''montage'' or ''video''. Current value: ''%s''.', visualizationFormat);
end

% Initialize varargout
varargout{1} = [];   % Impulse responses if (visualizeOuterSegmentFilters is true)
varargout{2} = [];   % noise-free responses
varargout{3} = [];   % the optical image optics
varargout{4} = [];   % the cone mosaic

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

%% Unique identifier for temporary block data prefix
blockDataPrefix = datestr(datetime('now'), 'yyyymmddTHHMMSS');

%% The constant params list
constantParamsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams,  rParams.backgroundParams, testDirectionParams};

colorModulationParamsNull = rParams.colorModulationParams;
colorModulationParamsNull.coneContrasts = [0 0 0]';
colorModulationParamsNull.contrast = 0;
    
%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
theProgram = mfilename;


%% Delete response instance files.
if (p.Results.delete)
    % Load the response and ancillary data
    paramsList = constantParamsList;
    paramsList{numel(paramsList)+1} = colorModulationParamsNull;
    fprintf('Removing null condition response instances.\n');
    rwObject.remove('responseInstances',paramsList,theProgram);

    paramsList = constantParamsList;
    paramsList{numel(paramsList)+1} = colorModulationParamsNull;   
    fprintf('Importing  ancillary data\n');
    ancillaryData = rwObject.read('ancillaryData',paramsList,theProgram);
    
    rParams = ancillaryData.rParams;
    parforConditionStructs = ancillaryData.parforConditionStructs;
    nParforConditions = length(parforConditionStructs); 
    for kk = nParforConditions:-1:1
            thisConditionStruct = parforConditionStructs{kk};
            colorModulationParamsTemp = rParams.colorModulationParams;
            colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
            colorModulationParamsTemp.contrast = thisConditionStruct.contrast;

            paramsList = constantParamsList;
            paramsList{numel(paramsList)+1} = colorModulationParamsTemp;    
            
            fprintf('Deleting condition %d/%d response instances\n', kk, nParforConditions);
            rwObject.remove('responseInstances',paramsList,theProgram);
    end % kk
end % if (p.Results.delete)

%% Maybe we are just visualizing the mosaic - not computing responses
if ((~p.Results.compute) && ((p.Results.visualizeMosaic) || (p.Results.computeMosaic)))
    if (p.Results.computeMosaic)
        % Create the cone mosaic
        theMosaic = colorDetectConeMosaicConstruct(rParams.mosaicParams, ...
            'visualizeMosaic', p.Results.visualizeMosaic);
        % Save cone mosaic
        coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
        rwObject.write('coneMosaic', theMosaic, coneParamsList, theProgram, 'type', 'mat');
    else
        % Load a previously saved cone mosaic
        fprintf('Loading a previously saved cone mosaic\n');
        coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
        theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
        theMosaic.displayInfo();
    end
    
    if (p.Results.visualizeMosaic) && (isa(theMosaic, 'coneMosaicHex'))
        hFig = theMosaic.visualizeGrid('generateNewFigure', true);
        coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
        data = 0;
        rwObject.write('coneMosaic', data, coneParamsList, theProgram, ...
            'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');   
    end
    
    if (p.Results.visualizeMosaic)
        % Return the mosaic
        varargout{4} = theMosaic;
    end
end

%% Maybe we are just visualizing the optics, not computing responses
if ((~p.Results.compute) && ((p.Results.visualizeOptics)) )

    % Load previously-generated optics
    fprintf('Loading previously generated optics\n');
    oiParamsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams};
    theOI = rwObject.read('opticalImage', oiParamsList, theProgram, 'type', 'mat');    
    
    % Return the OI (for visualization purposes)
    varargout{3} = theOI;
end

%% The computing happens here, if we are doing it
if (p.Results.compute)
    
    % Get the current time and date
    currentDate = datestr(datetime('now'), 'yyyymmddTHHMMSS');

    if (p.Results.computeOptics)
        % Create the optics
        if (ismember(rParams.oiParams.opticsModel, availableCustomWvfOpticsModels))
            [theOI, Zcoeffs] = colorDetectOpticalImageConstruct(rParams.oiParams);
        else
            theOI = colorDetectOpticalImageConstruct(rParams.oiParams);
        end

        % Save the optical image
        oiParamsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams};
        rwObject.write('opticalImage', theOI, oiParamsList, theProgram, 'type', 'mat');

        % Save a copy of the optical image with the current date
        % This can be useful if the optical image is overwritten 
        rwObject.write(sprintf('opticalImage_%s', currentDate), theOI, oiParamsList, theProgram, 'type', 'mat');

        % Save wvf Zcoeffs (if we are using custom wvf optics)
        if (ismember(rParams.oiParams.opticsModel, availableCustomWvfOpticsModels))
            rwObject.write('wvfZcoeffs', Zcoeffs, oiParamsList, theProgram, 'type', 'mat');
            rwObject.write(sprintf('wvfZcoeffs_%s', currentDate), Zcoeffs, oiParamsList, theProgram, 'type', 'mat');
        end
    else
        % Load previously-generated optics
        fprintf('Loading previously generated optics\n');
        oiParamsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams};
        theOI = rwObject.read('opticalImage', oiParamsList, theProgram, 'type', 'mat');    
    end
    
    % Return the OI (for visualization purposes)
    varargout{3} = theOI;

    if (p.Results.computeMosaic)
        % Create the cone mosaic
        theMosaic = colorDetectConeMosaicConstruct(rParams.mosaicParams, ...
            'visualizeMosaic', p.Results.visualizeMosaic);
        
        % Save cone mosaic
        coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
        rwObject.write('coneMosaic', theMosaic, coneParamsList, theProgram, 'type', 'mat');
        
        % Save a copy of the cone mosaic with the current date
        % This can be useful if the cone mosaic is overwritten 
        rwObject.write(sprintf('coneMosaic_%s', currentDate), theMosaic, coneParamsList, theProgram, 'type', 'mat');
    
    else
        % Load a previously saved cone mosaic
        fprintf('Loading a previously saved cone mosaic\n');
        coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
        theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
        theMosaic.displayInfo();
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
            
        case 'pedestalIncrements'
            if (strcmp(testDirectionParams.contrastScale, 'linear'))
                testContrasts = linspace(testDirectionParams.lowContrast, testDirectionParams.highContrast, testDirectionParams.nContrastsPerDirection);
            else
                testContrasts = logspace(log10(testDirectionParams.lowContrast), log10(testDirectionParams.highContrast), testDirectionParams.nContrastsPerDirection);
            end
            testConeContrasts = ones(3,1);
            
        otherwise
            error('Unknown test instance type passed: ''%s''.', testDirectionParams.instanceType);
    end
    
    parforConditionStructs = responseGenerationParforConditionStructsGenerate(testConeContrasts,testContrasts);
    nParforConditions = length(parforConditionStructs);
    nParforTrials = computeTrialBlocks(p.Results.ramPercentageEmployed, testDirectionParams.trialsNum, numel(theMosaic.pattern), numel(theMosaic.pattern(theMosaic.pattern>1)), rParams.temporalParams, theMosaic.integrationTime, p.Results.displayTrialBlockPartitionDiagnostics, p.Results.employStandardHostComputerResources);
    nParforTrialBlocks = numel(nParforTrials);
    parforRanSeeds = randi(1000000, nParforConditions, nParforTrialBlocks)+1;
    parforRanSeedsNoStim = randi(1000000,1, nParforTrialBlocks)+1;
    
    % Create parallel pool
    poolOBJ = gcp('nocreate');
    if isempty(poolOBJ)
        poolOBJ = parpool(parforWorkersNum);
    else
       delete(poolOBJ);
       poolOBJ = parpool(parforWorkersNum);
    end
    
    % Generate data for the no stimulus condition
    fprintf('<strong>[%02d] Computing the null stimulus response ... </strong>', 0);
    stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.2f', ...
        colorModulationParamsNull.coneContrasts(1), colorModulationParamsNull.coneContrasts(2), colorModulationParamsNull.coneContrasts(3), colorModulationParamsNull.contrast);
    
    paramsList = constantParamsList;
    paramsList{numel(paramsList)+1} = colorModulationParamsNull;
    
    
    % Start timing the null reponses computation
    tBegin = clock;
    
    % Parfor over trial blocks
    parfor (trialBlock = 1:nParforTrialBlocks, parforWorkersNum)
        % Get the parallel pool worker ID
        t = getCurrentTask();
        if (~isempty(p.Results.workerID))  
            workerID = t.ID;
        else
            workerID = [];
        end

        [tmpData{trialBlock}, osImpulseResponseFunctionsFromNullStimulus{trialBlock}, osImpulseReponseFunctionTimeAxis{trialBlock}, meanCurrentsFromNullStimulus{trialBlock}]  = ...
            colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nParforTrials(trialBlock), ...
                rParams.spatialParams, rParams.backgroundParams, colorModulationParamsNull, rParams.temporalParams, theOI, theMosaic, ...
                'centeredEMPaths', p.Results.centeredEMPaths, ...
                'osImpulseResponseFunctions', [], ...  % pass empty array, to compute the IR filters based on the null stimulus 
                'osMeanCurrents', [], ...              % pass empty array, to compute the steady-state current based on the null stimulus 
                'seed', parforRanSeedsNoStim(1, trialBlock), ...
                'workerID', workerID,...
                'displayTrialBlockPartitionDiagnostics', p.Results.displayTrialBlockPartitionDiagnostics, ...
                'computePhotocurrentResponseInstances', p.Results.computePhotocurrentResponseInstances, ...
                'computeNoiseFreeSignals', true, ....
                'visualizeSpatialScheme', (p.Results.visualizeSpatialScheme & (trialBlock == 1)), ...
                'visualizeOIsequence', false, ...
                'visualizeMosaicWithFirstEMpath', false, ...
                'paramsList', paramsList...
                );
                
            % save data temporarily
            rwObject.write(sprintf('%s_blockData_%d', blockDataPrefix, trialBlock), tmpData{trialBlock}, paramsList,theProgram);
            % Empty it to save space
            tmpData{trialBlock} = struct();
    end % parfor trialBlock
        
    % Form noStimData structure
    noStimData = struct(...
        'IBIOColorDetectSnapshot', p.Results.IBIOColorDetectSnapshot, ...
        'testContrast', colorModulationParamsNull.contrast, ...
        'testConeContrasts', colorModulationParamsNull.coneContrasts, ...
        'stimulusLabel', stimulusLabel, ...
        'osImpulseResponses', osImpulseResponseFunctionsFromNullStimulus{1},  ...
        'osImpulseResponseTimeAxis', osImpulseReponseFunctionTimeAxis{1}, ...
        'osMeanCurrents', meanCurrentsFromNullStimulus{1}, ...
        'responseInstanceArray', [], ...
        'noiseFreeIsomerizations', [], ...
        'noiseFreePhotocurrents', []);
    
    % Reassemble the responseInstanceArray for all trials from the blocked trials temporary datafiles
    [noStimData.responseInstanceArray, noStimData.noiseFreeIsomerizations, noStimData.noiseFreePhotocurrents] = ...
            formResponseInstanceArrayForAllTrials(rwObject, nParforTrialBlocks, paramsList, theProgram, blockDataPrefix, p.Results.displayTrialBlockPartitionDiagnostics);
    
    % Write the no cone contrast data and some extra facts we need
    rwObject.write('responseInstances',noStimData,paramsList,theProgram);
    
    % Clear to save RAM
    clear 'noStimData';
    
    % Save the other data we need for use by the classifier preprocessing subroutine
    ancillaryData = struct(...
        'testConeContrasts', testConeContrasts, ...
        'testContrasts', testContrasts, ...
        'rParams', rParams, ...
        'instanceParams', testDirectionParams);
    ancillaryData.parforConditionStructs = parforConditionStructs;
    rwObject.write('ancillaryData',ancillaryData,paramsList,theProgram);
    
    tEnd = clock;
    timeLapsed = etime(tEnd,tBegin);
    fprintf('<strong>Finished with the null stimulus response computations in %2.1f minutes. </strong>\n', timeLapsed/60);
    
    %% Generate data for all the examined stimuli
    % Loop over color directions
    stimDataForValidation = cell(nParforConditions,1);
    
    % Start timing the stim reponses computation
    tBegin = clock;
    
    for kk = nParforConditions:-1:1
        fprintf('<strong>[%02d] Computing responses for condition %d/%d ...  </strong> \n', kk, kk,nParforConditions);
        
        thisConditionStruct = parforConditionStructs{kk};
        colorModulationParamsTemp = rParams.colorModulationParams;
        colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
        colorModulationParamsTemp.contrast = thisConditionStruct.contrast;
        
        paramsList = constantParamsList;
        paramsList{numel(paramsList)+1} = colorModulationParamsTemp;
            
        % Make noisy instances for each contrast
        stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.5f',...
            colorModulationParamsTemp.coneContrasts(1), colorModulationParamsTemp.coneContrasts(2), colorModulationParamsTemp.coneContrasts(3), colorModulationParamsTemp.contrast);
    
        osImpulseResponses = {};
        osMeanCurrents = {};
        
        % Parfor over blocks of trials
        parfor (trialBlock = 1:nParforTrialBlocks, parforWorkersNum) 
            % Get the parallel pool worker ID
            t = getCurrentTask();
            
            if (~isempty(p.Results.workerID))  
                workerID = t.ID;
            else
                workerID = [];
            end
            
            [tmpData{trialBlock}, osImpulseResponses{trialBlock}, osImpulseReponseFunctionTimeAxis{trialBlock}, osMeanCurrents{trialBlock}] = ...
                colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nParforTrials(trialBlock), ...
                    rParams.spatialParams, rParams.backgroundParams, colorModulationParamsTemp, rParams.temporalParams, theOI, theMosaic, ...
                    'centeredEMPaths', p.Results.centeredEMPaths, ...
                    'osImpulseResponseFunctions', osImpulseResponseFunctionsFromNullStimulus{1}, ...  % use the IR filters computed from the null stimulus
                    'osMeanCurrents', meanCurrentsFromNullStimulus{1}, ...                            % use the to steady-state currents computed from the null stimulus  
                    'seed', parforRanSeeds(kk,trialBlock), ...
                    'workerID', workerID, ...
                    'displayTrialBlockPartitionDiagnostics', p.Results.displayTrialBlockPartitionDiagnostics, ...
                    'computePhotocurrentResponseInstances', p.Results.computePhotocurrentResponseInstances, ...
                    'computeNoiseFreeSignals', true, ...
                    'visualizeSpatialScheme', (p.Results.visualizeSpatialScheme & (trialBlock == 1)), ...
                    'visualizeOIsequence', (p.Results.visualizeOIsequence & (trialBlock == 1) & (kk == nParforConditions)), ...
                    'visualizeMosaicWithFirstEMpath', (p.Results.visualizeMosaicWithFirstEMpath & (trialBlock == 1)), ...
                    'paramsList', paramsList);
            
            % Save data temporarily
            rwObject.write(sprintf('%s_blockData_%d', blockDataPrefix, trialBlock), tmpData{trialBlock}, paramsList,theProgram);
            % Empty it to save space
            tmpData{trialBlock} = struct();
        end % parfor trialBlock
        
        % Form stimData structure
        stimData = struct(...
            'IBIOColorDetectSnapshot', p.Results.IBIOColorDetectSnapshot, ...
            'testContrast', colorModulationParamsTemp.contrast, ...
            'testConeContrasts', colorModulationParamsTemp.coneContrasts, ...
            'stimulusLabel', stimulusLabel, ...
            'osImpulseResponses', osImpulseResponses{1}, ...
            'osImpulseResponseTimeAxis', osImpulseReponseFunctionTimeAxis{1}, ...
            'osMeanCurrents', osMeanCurrents{1}, ...
            'responseInstanceArray', [], ...
            'noiseFreeIsomerizations', [], ...
            'noiseFreePhotocurrents', []);
        
        % Reassemble the responseInstanceArray for all trials from the blocked trials temporary datafiles
        [stimData.responseInstanceArray, stimData.noiseFreeIsomerizations, stimData.noiseFreePhotocurrents] = ...
            formResponseInstanceArrayForAllTrials(rwObject, nParforTrialBlocks, paramsList, theProgram, blockDataPrefix, p.Results.displayTrialBlockPartitionDiagnostics);

        % Save some data for validation in first loop
        if (kk == 1)
            s = struct();
            savedTrial = 1;
            s.theMosaicIsomerizations = squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(savedTrial,:,:,:));
            if (isempty(stimData.responseInstanceArray.theMosaicPhotocurrents))
                s.theMosaicPhotocurrents = [];
            else
                s.theMosaicPhotocurrents = squeeze(stimData.responseInstanceArray.theMosaicPhotocurrents(savedTrial,:,:,:));
            end
            s.theMosaicEyeMovements = squeeze(stimData.responseInstanceArray.theMosaicEyeMovements(savedTrial,:,:));
            s.timeAxis = stimData.responseInstanceArray.timeAxis;
            s.photocurrentTimeAxis = stimData.responseInstanceArray.photocurrentTimeAxis;
            stimDataForValidation = s;
        end
        
        % Save data for this color direction/contrast pair
        rwObject.write('responseInstances',stimData,paramsList,theProgram);
        
        % Clear to save RAM
        clear 'stimData';
    end % for kk = conditions
    
    tEnd = clock;
    timeLapsed = etime(tEnd,tBegin);
    fprintf('<strong>Finished response computation for all conditions in %2.3f hours. </strong>\n', timeLapsed/60/60);
    
    %% Validation data
    %
    % Since we want to store some stimData and that isn't saved outside the
    % parfor loop (and is in fact a bit hard to save outside the parfor
    % loop)
    if (nargout > 0)
        % Reload the no-stim data
        paramsList = constantParamsList;
        paramsList{numel(paramsList)+1} = colorModulationParamsNull;   
        noStimData = rwObject.read('responseInstances',paramsList,theProgram);
    
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
        clear 'noStimData';
        
        validationData.noStimData = noStimValidationData;
        validationData.stimData = stimDataForValidation;
        extraData.ancillaryData = ancillaryData;
        extraData.p.Results = p.Results;
    end
end

%% Visualize
if ((p.Results.visualizeResponses || p.Results.visualizeOuterSegmentFilters))

    if (p.Results.visualizeResponses)
        % Load the mosaic
        coneParamsList = {rParams.topLevelDirParams, rParams.mosaicParams};
        theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');

        % Load the response and ancillary data
        paramsList = constantParamsList;
        paramsList{numel(paramsList)+1} = colorModulationParamsNull;   
        fprintf('Importing NOSTIM data and ancillary data\n');
        noStimData = rwObject.read('responseInstances',paramsList,theProgram);
        ancillaryData = rwObject.read('ancillaryData',paramsList,theProgram);
    
         % Save noise-free data (to be returned to the user)
        noiseFreeResponses.noStimData.noiseFreeIsomerizations = noStimData.noiseFreeIsomerizations;
        noiseFreeResponses.noStimData.noiseFreePhotocurrents  = noStimData.noiseFreePhotocurrents;
    
        % How many istances to visualize
        instancesToVisualize = size(noStimData.responseInstanceArray.theMosaicIsomerizations,1); % 5;

        % Only keep the data we will visualize
        if (instancesToVisualize ~= size(noStimData.responseInstanceArray.theMosaicIsomerizations,1))
            idx = min([instancesToVisualize size(noStimData.responseInstanceArray.theMosaicIsomerizations,1)]);
            noStimData.responseInstanceArray.theMosaicIsomerizations = noStimData.responseInstanceArray.theMosaicIsomerizations(1:idx,:,:);
            if (~isempty(noStimData.responseInstanceArray.theMosaicPhotocurrents))
                noStimData.responseInstanceArray.theMosaicPhotocurrents  = noStimData.responseInstanceArray.theMosaicPhotocurrents(1:idx,:,:);
            end
            noStimData.responseInstanceArray.theMosaicEyeMovements = noStimData.responseInstanceArray.theMosaicEyeMovements(1:idx,:,:);
        end
        
        rParams = ancillaryData.rParams;
        parforConditionStructs = ancillaryData.parforConditionStructs;
        nParforConditions = length(parforConditionStructs); 
        for kk = nParforConditions:-1:1
            if (~isempty(p.Results.visualizedConditionIndices)) && (~ismember(kk, p.Results.visualizedConditionIndices))
                fprintf('Skipping condition %d/%d STIM data for visualization\n', kk, nParforConditions);
                continue;
            end
            fprintf('Importing condition %d/%d STIM data for visualization\n', kk, nParforConditions);
            thisConditionStruct = parforConditionStructs{kk};
            colorModulationParamsTemp = rParams.colorModulationParams;
            colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
            colorModulationParamsTemp.contrast = thisConditionStruct.contrast;

            paramsList = constantParamsList;
            paramsList{numel(paramsList)+1} = colorModulationParamsTemp;    
            stimData = rwObject.read('responseInstances',paramsList,theProgram);
            
            % Save noise-free data (to be returned to the user)
            noiseFreeResponses.stimData{kk}.noiseFreeIsomerizations = stimData.noiseFreeIsomerizations;
            noiseFreeResponses.stimData{kk}.noiseFreePhotocurrents  = stimData.noiseFreePhotocurrents;
        
            if (instancesToVisualize ~= size(noStimData.responseInstanceArray.theMosaicIsomerizations,1))
                % Only keep the data we will visualize
                stimData.responseInstanceArray.theMosaicIsomerizations = stimData.responseInstanceArray.theMosaicIsomerizations(1:idx,:,:);
                if (~isempty(stimData.responseInstanceArray.theMosaicPhotocurrents))
                    stimData.responseInstanceArray.theMosaicPhotocurrents = stimData.responseInstanceArray.theMosaicPhotocurrents(1:idx,:,:);
                end
                stimData.responseInstanceArray.theMosaicEyeMovements = stimData.responseInstanceArray.theMosaicEyeMovements(1:idx,:,:);
            end
            
            % Visualization routines work only for coneMosaicHex
            if (p.Results.generatePlots) && (isa(theMosaic, 'coneMosaicHex'))
                % Call the figures generation routine
                hFigsInfo = visualizeResponseInstances(theMosaic, stimData, noStimData, p.Results.visualizeOuterSegmentFilters, p.Results.visualizedResponseNormalization, kk, nParforConditions, p.Results.visualizationFormat);

                % Save figures produced
                theProgram = mfilename;
                rwObject = IBIOColorDetectReadWriteBasic;
                data = 0;

                for k = 1:numel(hFigsInfo)
                    if ~isempty(hFigsInfo{k}.hFig)
                        fileName = hFigsInfo{k}.filename;
                        rwObject.write(fileName, data, paramsList, theProgram, ...
                            'type', 'NicePlotExportPDF', 'FigureHandle', hFigsInfo{k}.hFig, 'FigureType', 'pdf');
                    end
                end
            end
        end
        
        varargout{2} = noiseFreeResponses;
    end
    
    if (p.Results.visualizeOuterSegmentFilters)
        % Load the response data
        paramsList = constantParamsList;
        paramsList{numel(paramsList)+1} = colorModulationParamsNull;   
        fprintf('Importing NOSTIM data\n');
        ancillaryData = rwObject.read('ancillaryData',paramsList,theProgram);
        rParams = ancillaryData.rParams;
        parforConditionStructs = ancillaryData.parforConditionStructs;
        nParforConditions = length(parforConditionStructs);   
        noStimData = rwObject.read('responseInstances',paramsList,theProgram);
        noStimImpulseResponses = noStimData.osImpulseResponses;
        osImpulseResponseTimeAxis = noStimData.osImpulseResponseTimeAxis;
        
        varargout{1} = noStimImpulseResponses;
        for kk = 1:nParforConditions
            fprintf('Importing STIM data for condition %d/%d\n', kk, nParforConditions);
            thisConditionStruct = parforConditionStructs{kk};
            colorModulationParamsTemp = rParams.colorModulationParams;
            colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
            colorModulationParamsTemp.contrast = thisConditionStruct.contrast;

            paramsList = constantParamsList;
            paramsList{numel(paramsList)+1} = colorModulationParamsTemp;    
            stimData = rwObject.read('responseInstances',paramsList,theProgram);
            stimImpulseResponses(kk,:,:) = stimData.osImpulseResponses;
        end % kk
        
        visualizeOuterSegmentImpulseResponseFunctions(osImpulseResponseTimeAxis, noStimImpulseResponses, stimImpulseResponses, constantParamsList);
    end
end

end

function [responseInstanceArray, noiseFreeIsomerizations, noiseFreePhotocurrents] = formResponseInstanceArrayForAllTrials(rwObject, nParforTrialBlocks, paramsList, theProgram, blockDataPrefix, displayTrialBlockPartitionDiagnostics)
    if (displayTrialBlockPartitionDiagnostics)
        fprintf('Importing data from %d blocks to build up the responseInstanceArray.\n', nParforTrialBlocks);
    end
    % Load all the block data
    for trialBlock = 1:nParforTrialBlocks
        blockDataFileName = sprintf('%s_blockData_%d', blockDataPrefix, trialBlock);
        tmpData = rwObject.read(blockDataFileName, paramsList, theProgram, 'type', 'mat');
        if (trialBlock == 1)
            noiseFreeIsomerizations = tmpData.noiseFreeIsomerizations;
            noiseFreePhotocurrents = tmpData.noiseFreePhotocurrents;
            responseInstanceArray = struct();
            responseInstanceArray.timeAxis = tmpData.responseInstanceArray.timeAxis;
            responseInstanceArray.photocurrentTimeAxis = tmpData.responseInstanceArray.photocurrentTimeAxis;
            responseInstanceArray.theMosaicEyeMovements = [];
            responseInstanceArray.theMosaicIsomerizations = [];
            responseInstanceArray.theMosaicPhotocurrents = [];
        end
        responseInstanceArray.theMosaicEyeMovements   = cat(1, responseInstanceArray.theMosaicEyeMovements,   tmpData.responseInstanceArray.theMosaicEyeMovements);
        responseInstanceArray.theMosaicIsomerizations = cat(1, responseInstanceArray.theMosaicIsomerizations, tmpData.responseInstanceArray.theMosaicIsomerizations);
        responseInstanceArray.theMosaicPhotocurrents  = cat(1, responseInstanceArray.theMosaicPhotocurrents,  tmpData.responseInstanceArray.theMosaicPhotocurrents);
        % Delete the block data file
        rwObject.remove(blockDataFileName,paramsList,theProgram);
    end % trialBlock
    if (displayTrialBlockPartitionDiagnostics)
        fprintf('Done with building up the responseInstanceArray\n');
    end
end

function nParforTrials = computeTrialBlocks(ramPercentageEmployed, nTrials, coneMosaicPatternSize, coneMosaicActivePatternSize, temporalParams, integrationTime, displayTrialBlockPartitionDiagnostics, employStandardHostComputerResources)      
    % Determine system resources
    [numberOfWorkers, ramSizeGBytes, sizeOfDoubleInBytes] = determineSystemResources(employStandardHostComputerResources);
   
    if (numberOfWorkers == 1)
        nParforTrials = nTrials;
        return;
    end
    
    % Ensure ramPercentageEmployed is in [0.05 1]
    ramPercentageEmployed = max([0.05 min([1 ramPercentageEmployed])]);
    
    ramSizeGBytes = ramPercentageEmployed * ramSizeGBytes;
    
    % Subtract RAM used by the OS
    ramUsedByOSGBytes = 1.2;
    ramSizeGBytesAvailable = ramSizeGBytes - ramUsedByOSGBytes;
    
    % Compute sizes of the large players
    if (isnan(temporalParams.windowTauInSeconds))
        [stimulusTimeAxis, ~, ~] = squareTemporalWindowCreate(temporalParams);
    else
        [stimulusTimeAxis, ~, ~] = gaussianTemporalWindowCreate(temporalParams);
    end

    
    if (numel(stimulusTimeAxis) == 1)
        stimulusSamplingInterval  = integrationTime;
    else
        stimulusSamplingInterval = stimulusTimeAxis(2)-stimulusTimeAxis(1);
    end
    
    eyeMovementsNumPerOpticalImage = stimulusSamplingInterval/integrationTime;
    emPathLength = round(eyeMovementsNumPerOpticalImage*numel(stimulusTimeAxis));
    
    % estimate sizes of the various matrices used
    trialBlockSize = floor(nTrials/numberOfWorkers);
    totalRAM = ramSizeGBytesAvailable;
    allowedRAMcompression = 0.9;
    while (totalRAM > allowedRAMcompression*ramSizeGBytesAvailable)
        totalMemoryPerWorker = 2*(coneMosaicPatternSize*trialBlockSize+coneMosaicActivePatternSize*emPathLength*trialBlockSize)*sizeOfDoubleInBytes/(1024^3);
        totalRAM = totalMemoryPerWorker * numberOfWorkers;
        trialBlockSize = trialBlockSize-1;
    end
    
    trialBlockSize = min([nTrials max([1 trialBlockSize])]);
    nParforTrialBlocks = ceil(nTrials / trialBlockSize);
    
    if (nParforTrialBlocks < numberOfWorkers)
        nParforTrialBlocks = numberOfWorkers;
        trialBlockSize = max([1 floor(nTrials/nParforTrialBlocks)]);
        totalMemoryPerWorker = 2*(coneMosaicPatternSize*trialBlockSize+coneMosaicActivePatternSize*emPathLength*trialBlockSize)*sizeOfDoubleInBytes/(1024^3);
    end
    
    totalMemoryUsed = numberOfWorkers * totalMemoryPerWorker;
    if (totalMemoryUsed < 0.01*ramSizeGBytesAvailable/numberOfWorkers)
         % Just use one processor - faster most of the times
         nParforTrials(1) = nTrials;
         nParforTrialBlocks = 1;
         trialBlockSize = nTrials;
    elseif (totalMemoryUsed < ramSizeGBytesAvailable/numberOfWorkers)
          nParforTrialBlocks = numberOfWorkers;
          trialBlockSize = floor(nTrials/nParforTrialBlocks);
          for kk = 1:nParforTrialBlocks-1
              nParforTrials(kk) = trialBlockSize;
          end
          nParforTrials(nParforTrialBlocks) = nTrials - trialBlockSize*nParforTrialBlocks;
          totalMemoryPerWorker = 2*(coneMosaicPatternSize*trialBlockSize+coneMosaicActivePatternSize*emPathLength*trialBlockSize)*sizeOfDoubleInBytes/(1024^3);
          totalMemoryUsed = totalMemoryPerWorker * numberOfWorkers;
    end

     
    
    nParforTrialsTmp = ones(1, nParforTrialBlocks) * trialBlockSize;
    remainingTrials = nTrials - sum(nParforTrialsTmp);
    if (remainingTrials > 0) 
        if remainingTrials > round(trialBlockSize*0.2)
            nParforTrialsTmp(end+1) = remainingTrials;
        else
            nParforTrialsTmp(end) = nParforTrialsTmp(end) + remainingTrials;
        end
    end
    
    accumTrials = 0;
    k = 1; keepGoing = true;
    while (k <= numel(nParforTrialsTmp)) && (keepGoing)
        if (accumTrials+nParforTrialsTmp(k) > nTrials)
            if (nTrials > accumTrials)
                nParforTrials(k) = nTrials - accumTrials;
            end
            keepGoing = false;
        else
            nParforTrials(k) = nParforTrialsTmp(k);
            accumTrials = accumTrials+nParforTrials(k);
            k = k + 1;
        end
    end

    if (sum(nParforTrials) ~= nTrials)
        nParforTrials
        nTrials
        trialBlockSize
        nParforTrialBlocks
        error('Error in logic of trial partitioning.')
    end
    
    fprintf('<strong> %d workers, system RAM = %2.1fGBytes </strong> \n', numberOfWorkers, ramSizeGBytes), ...
    fprintf('<strong> %d trials partitioned in %d blocks, each with %d trials (last has %d trials) </strong>\n', nTrials, numel(nParforTrials), nParforTrials(1), nParforTrials(end));
    fprintf('<strong> RAM used : %2.1f GBytes (per worker), total: %2.1f GBytes </strong> \n\n', totalMemoryPerWorker, totalMemoryUsed);
    
end

