function [responseInstanceArray,noiseFreeIsomerizations,noiseFreePhotocurrents, osImpulseResponseFunctions, osMeanCurrents] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic, varargin)
% [responseInstanceArray,noiseFreeIsomerizations,noiseFreePhotocurrents] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic)
% 
% Construct an array of nTrials response instances given the
% simulationTimeStep, spatialParams, temporalParams, theOI, theMosaic.
%
% The noise free isomerizations response is returned for the first
% instance. NOTE THAT BECAUSE EYE MOVEMENT PATHS DIFFER ACROSS INSTANCES,
% THIS IS NOT THE SAME FOR EVERY INSTANCE.]
%
% Key/value pairs
%  'seed' - value (default 1). Random number generator seed
%   'workerID' - (default empty).  If this field is non-empty, the progress of
%            the computation is printed in the command window along with the
%            workerID (from a parfor loop).
%   'trialBlockSize' - How many trials each trialBlock should have. Default: [], which results in nTrials (no blocking). 
%               This only has an effect with @coneMosaicHex mosaics and when nTrials>1 and it is useful with 
%               large mosaics x lots of trials, in which case the absorptions matrix does not fit in the RAM.
%               If set to -1, the number of trial blocks is computed automatically based on the number of cores and system RAM.
%  'centeredEMPaths' - true/false (default false) 
%               Controls wether the eye movement paths start at (0,0) (default) or wether they are centered around (0,0)
%  'osImpulseResponseFunctions' - the LMS impulse response filters to be used (default: [])
%  'osMeanCurrents' - the steady-state LMS currents caused by the mean absorption LMS rates
%  'computeNoiseFreeSignals' - true/false (default true) wether to compute the noise free isomerizations and photocurrents
%  'useSinglePrecision' - true/false (default true) use single precision to represent isomerizations and photocurrent
%  'displayTrialBlockPartitionDiagnostics' - true/false (default false)
%               Controls wether to display how trials are partitioned into trial blocks
% 7/10/16  npc Wrote it.

%% Parse arguments
p = inputParser;
p.addParameter('seed',1, @isnumeric);   
p.addParameter('workerID', [], @isnumeric);
p.addParameter('useSinglePrecision',true,@islogical);
p.addParameter('trialBlockSize', [], @isnumeric);
p.addParameter('centeredEMPaths',false, @islogical); 
p.addParameter('osImpulseResponseFunctions', [], @isnumeric);
p.addParameter('osMeanCurrents', [], @isnumeric);
p.addParameter('computeNoiseFreeSignals', true, @islogical);
p.addParameter('displayTrialBlockPartitionDiagnostics', false, @islogical);

p.parse(varargin{:});
currentSeed = p.Results.seed;
trialBlockSize = p.Results.trialBlockSize;
displayTrialBlockPartitionDiagnostics = p.Results.displayTrialBlockPartitionDiagnostics;

%% Start computation time measurement
tic

%% Get pupil size out of OI, which is sometimes needed by colorSceneCreate
oiParamsTemp.pupilDiamMm = 1000*opticsGet(oiGet(theOI,'optics'),'aperture diameter');

%% Create the background scene
theBaseColorModulationParams = colorModulationParams;
theBaseColorModulationParams.coneContrasts = [0 0 0]';
theBaseColorModulationParams.contrast = 0;
backgroundScene = colorSceneCreate(spatialParams, backgroundParams, ...
    theBaseColorModulationParams, oiParamsTemp);

%% Create the modulated scene
modulatedScene = colorSceneCreate(spatialParams, backgroundParams, ...
    colorModulationParams, oiParamsTemp);

%% Compute the background OI
oiBackground = theOI;
oiBackground = oiCompute(oiBackground, backgroundScene);

%% Compute the modulated OI
oiModulated = theOI;
oiModulated = oiCompute(oiModulated, modulatedScene);

%% Compute size of optical image (used in determining trialBlockSize later on)
s = whos('oiModulated');
oiSizeBytes = s.bytes;

%% Generate the stimulus modulation function
[stimulusTimeAxis, stimulusModulationFunction, ~] = gaussianTemporalWindowCreate(temporalParams);

%% Compute the oiSequence
theOIsequence = oiSequence(oiBackground, oiModulated, stimulusTimeAxis, ...
                                stimulusModulationFunction, 'composition', 'blend');

% Clear oiModulated and oiBackground to save space
varsToClear = {'oiBackground', 'oiModulated'};
clear(varsToClear{:});

%% Generate eye movement paths for all instances
eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theMosaic.integrationTime);
theEMpaths = colorDetectMultiTrialEMPathGenerate(...
    theMosaic, nTrials, eyeMovementsNum, temporalParams.emPathType, ...
    'centeredEMPaths', p.Results.centeredEMPaths, ...
    'seed', currentSeed);

% Determine optimal trialBlocks
computeTrialBlockSizeForParforLoop(nTrials, numel(theMosaic.pattern), numel(theMosaic.pattern(theMosaic.pattern>1)), oiSizeBytes, eyeMovementsNum, displayTrialBlockPartitionDiagnostics);
if (trialBlockSize == -1)
    trialBlockSize = computeTrialBlockSizeForParforLoop(nTrials, numel(theMosaic.pattern), numel(theMosaic.pattern(theMosaic.pattern>1)), oiSizeBytes, eyeMovementsNum, displayTrialBlockPartitionDiagnostics);
end

[isomerizations, photocurrents, osImpulseResponseFunctions, osMeanCurrents] = ...
     theMosaic.computeForOISequence(theOIsequence, ...
                    'seed', currentSeed, ...
                    'emPaths', theEMpaths, ...
                    'interpFilters', p.Results.osImpulseResponseFunctions, ...
                    'meanCur', p.Results.osMeanCurrents, ...
                    'trialBlockSize', trialBlockSize, ...
                    'currentFlag', true, ...
                    'workDescription', 'noisy responses', ...
                    'workerID', p.Results.workerID);
isomerizationsTimeAxis = theMosaic.timeAxis + theOIsequence.timeAxis(1);
photoCurrentTimeAxis = isomerizationsTimeAxis;
 
%% Remove unwanted portions of the responses
isomerizationsTimeIndicesToKeep = find(abs(isomerizationsTimeAxis-temporalParams.secondsToIncludeOffset) <= temporalParams.secondsToInclude/2);
photocurrentsTimeIndicesToKeep = find(abs(photoCurrentTimeAxis-temporalParams.secondsToIncludeOffset) <= temporalParams.secondsToInclude/2);
if (isa(theMosaic , 'coneMosaicHex'))
     isomerizations = isomerizations(:,:,isomerizationsTimeIndicesToKeep);
     photocurrents = photocurrents(:,:,photocurrentsTimeIndicesToKeep);
else
     isomerizations = isomerizations(:,:,:,isomerizationsTimeIndicesToKeep);
     photocurrents = photocurrents(:,:,:,photocurrentsTimeIndicesToKeep);
end

%% Form the responseInstanceArray struct
if (p.Results.useSinglePrecision)
    responseInstanceArray.theMosaicIsomerizations = single(isomerizations);
    responseInstanceArray.theMosaicPhotocurrents = single(photocurrents);
else
    responseInstanceArray.theMosaicIsomerizations = isomerizations;
    responseInstanceArray.theMosaicPhotocurrents = photocurrents;
end
responseInstanceArray.theMosaicEyeMovements = theEMpaths(:,isomerizationsTimeIndicesToKeep,:);
responseInstanceArray.timeAxis = isomerizationsTimeAxis(isomerizationsTimeIndicesToKeep);
responseInstanceArray.photocurrentTimeAxis = photoCurrentTimeAxis(photocurrentsTimeIndicesToKeep);


if (~p.Results.computeNoiseFreeSignals)
    noiseFreeIsomerizations = [];
    noiseFreePhotocurrents = [];
    return;
end

%% Compute the noise-free isomerizations & photocurrents
if strcmp(temporalParams.emPathType, 'frozen')
    % use the first emPath if we have a frozen path
    theEMpaths = theEMpaths(1,:,:);
else
    % use all zeros path otherwise
    theEMpaths = theEMpaths(1,:,:)*0;
end

% Save original noise flags
originalIsomerizationNoiseFlag = theMosaic.noiseFlag;
originalPootocurrentNoiseFlag = theMosaic.os.noiseFlag;

% Set noiseFlags to none
theMosaic.noiseFlag = 'none';
theMosaic.os.noiseFlag = 'none';

% Compute the noise-free responses
[noiseFreeIsomerizations, noiseFreePhotocurrents] = ...
     theMosaic.computeForOISequence(theOIsequence, ...
                    'emPaths',theEMpaths, ...
                    'interpFilters', p.Results.osImpulseResponseFunctions, ...
                    'meanCur', p.Results.osMeanCurrents, ...
                    'currentFlag', true, ...
                    'workDescription', 'noise-free responses', ...
                    'workerID', p.Results.workerID);

%% Remove unwanted portions of the noise-free responses
if (isa(theMosaic , 'coneMosaicHex'))
     noiseFreeIsomerizations = noiseFreeIsomerizations(:,:,isomerizationsTimeIndicesToKeep);
     noiseFreePhotocurrents = noiseFreePhotocurrents(:,:,photocurrentsTimeIndicesToKeep);
else
     noiseFreeIsomerizations = noiseFreeIsomerizations(:,:,:,isomerizationsTimeIndicesToKeep);
     noiseFreePhotocurrents = noiseFreePhotocurrents(:,:,:,photocurrentsTimeIndicesToKeep);
end

% Restore noiseFlags to none
theMosaic.noiseFlag = originalIsomerizationNoiseFlag;
theMosaic.os.noiseFlag = originalPootocurrentNoiseFlag;

% Store
noiseFreeIsomerizations = squeeze(noiseFreeIsomerizations);
if (~isempty(noiseFreePhotocurrents))
    noiseFreePhotocurrents = squeeze(noiseFreePhotocurrents);
end

%% Report time taken
fprintf('Response instance array generation (%d instances) took %2.3f minutes to compute.\n', nTrials, toc/60);

end



function trialBlockSize = computeTrialBlockSizeForParforLoop(nTrials, coneMosaicPatternSize, coneMosaicActivePatternSize, oiSizeBytes, emPathLength, displayDiagnostics)
        
    % Determine system resources
    [numberOfWorkers, ramSizeGBytes, sizeOfDoubleInBytes] = determineSystemResources();

    % Subtract RAM used by the OS
    ramUsedByOSGBytes = 1.2;
    ramSizeGBytesAvailable = ramSizeGBytes - ramUsedByOSGBytes;
    
    % Compute sizes of the large players
    % mosaic.compute temporary products
    computeTempProductsSizeGBytes = coneMosaicPatternSize*sizeOfDoubleInBytes/(1024^3);
    
    absorptions_memsizeGBytes = nTrials*coneMosaicActivePatternSize*emPathLength*sizeOfDoubleInBytes/2/(1024^3);
    photocurrents_memsizeGBytes = absorptions_memsizeGBytes;
    
    obj_currents_memsizeGBytes = coneMosaicActivePatternSize*emPathLength*sizeOfDoubleInBytes/(1024^3);
    obj_absorptions_memsizeGBytes = coneMosaicActivePatternSize*sizeOfDoubleInBytes/(1024^3);
    
    % OIsizes (oiModulated, oiFixed, + mixture + currentOI + previousOI)
    oi_memsizeGBytes = (2 * oiSizeBytes)/(1024^3);
    
    % Compute trialBlocksSize
    totalMemoryPerWorker = max(...
        [ computeTempProductsSizeGBytes+obj_absorptions_memsizeGBytes+oi_memsizeGBytes  ...
          absorptions_memsizeGBytes+photocurrents_memsizeGBytes+obj_currents_memsizeGBytes...
        ]);
    singleTrialMemoryGBytes = numberOfWorkers * totalMemoryPerWorker;
    
    allowedRAMcompression = 1.0;
    trialBlockSize = round(allowedRAMcompression*ramSizeGBytesAvailable/singleTrialMemoryGBytes);
    if (trialBlockSize > nTrials)
        trialBlockSize = nTrials;
    end
    blockedTrialIndices = computeBlockedTrialIndices(trialBlockSize, nTrials);

    % Display nTrial partitioning into blocks
    if (displayDiagnostics)
        for iTrialBlock = 1:numel(blockedTrialIndices)
            trialIndicesForThisBlock = blockedTrialIndices{iTrialBlock};
            firstTrial = trialIndicesForThisBlock(1);
            lastTrial = trialIndicesForThisBlock(end);
            if (iTrialBlock == 1)
                firstTrialBlockSize = lastTrial-firstTrial+1;
            end
            if (iTrialBlock == numel(blockedTrialIndices))
                lastTrialBlockSize = lastTrial-firstTrial+1;
            end
            %fprintf('trialBlock%d contains %d trials: %04d - %04d\n', iTrialBlock, numel(trialIndicesForThisBlock), firstTrial, lastTrial);
        end
        fprintf('----------------------------------------------------------------------------------\n');
        fprintf('<strong>CoresNum = %d; SystemRAM = %2.2fGB; estimated peak RAM = %2.2fGB</strong>\n', numberOfWorkers, ramSizeGBytes, trialBlockSize*singleTrialMemoryGBytes+ramUsedByOSGBytes);
        fprintf('<strong>%d trials prganized in %d blocks, blockSize(1/last) = %d/%d</strong>\n', nTrials, numel(blockedTrialIndices), firstTrialBlockSize, lastTrialBlockSize);
        fprintf('----------------------------------------------------------------------------------\n');
    end
end
