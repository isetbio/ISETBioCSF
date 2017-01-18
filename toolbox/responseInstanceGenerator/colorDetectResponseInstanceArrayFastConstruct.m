function [responseInstanceArray,noiseFreeIsomerizations,noiseFreePhotocurrents] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic, varargin)
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
%  'useSinglePrecision' - true/false (default true) use single precision to represent isomerizations and photocurrent

% 7/10/16  npc Wrote it.

%% Parse arguments
p = inputParser;
p.addParameter('seed',1, @isnumeric);   
p.addParameter('workerID', [], @isnumeric);
p.addParameter('useSinglePrecision',true,@islogical);
p.addParameter('trialBlockSize', [], @isnumeric);
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

%% Generate the stimulus modulation function
[stimulusTimeAxis, stimulusModulationFunction, ~] = gaussianTemporalWindowCreate(temporalParams);

%% Compute the oiSequence
theOIsequence = oiSequence(oiBackground, oiModulated, stimulusTimeAxis, ...
                                stimulusModulationFunction, 'composition', 'blend');

%% Generate eye movement paths for all instances
eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theMosaic.integrationTime);
theEMpaths = colorDetectMultiTrialEMPathGenerate(theMosaic, nTrials, eyeMovementsNum, temporalParams.emPathType, 'seed', currentSeed);


% Determine optimal trialBlocks
computeTrialBlockSizeForParforLoop(nTrials, numel(theMosaic.pattern), numel(theMosaic.pattern(theMosaic.pattern>1)), eyeMovementsNum, displayTrialBlockPartitionDiagnostics);
if (trialBlockSize == -1)
    trialBlockSize = computeTrialBlockSizeForParforLoop(nTrials, numel(theMosaic.pattern), numel(theMosaic.pattern(theMosaic.pattern>1)), eyeMovementsNum, displayTrialBlockPartitionDiagnostics);
end

[isomerizations, photocurrents] = ...
     theMosaic.computeForOISequence(theOIsequence, ...
                    'seed', currentSeed, ...
                    'emPaths', theEMpaths, ...
                    'trialBlockSize', trialBlockSize, ...
                    'currentFlag', true, ...
                    'workDescription', 'noisy responses', ...
                    'workerID', p.Results.workerID);
isomerizationsTimeAxis = theMosaic.timeAxis + theOIsequence.timeAxis(1);
photoCurrentTimeAxis = isomerizationsTimeAxis;
 
%% Remove unwanted portions of the responses
isomerizationsTimeIndicesToKeep = find(abs(isomerizationsTimeAxis-temporalParams.secondsToIncludeOffset) <= temporalParams.secondsToInclude/2);
photocurrentsTimeIndicesToKeep = find(abs(photoCurrentTimeAxis-temporalParams.secondsToIncludeOffset) <= temporalParams.secondsToInclude/2);

%% Reshape if we're using a hex mosaic
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

%% Get noise-free responses
%
% Compute the noise-free isomerizations & photocurrents using the first emPath
theEMpaths = theEMpaths(1,:,:);

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
                    'currentFlag', true, ...
                    'workDescription', 'noise-free responses', ...
                    'workerID', p.Results.workerID);
                
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



function trialBlockSize = computeTrialBlockSizeForParforLoop(nTrials, coneMosaicPatternSize, coneMosaicActivePatternSize, emPathLength, displayDiagnostics)
        
    % Determine system resources
    [numberOfCores, ramSizeGBytes, sizeOfDoubleInBytes] = determineSystemResources();

    % Subtract RAM used by the OS
    ramUsedByOSGBytes = 3.2;
    ramSizeGBytesAvailable = ramSizeGBytes - ramUsedByOSGBytes;
    
    % Compute trialBlocksForParforLoop
    allowedRAMcompression = 0.8;
    singleTrialMemoryGBytes = 2*numberOfCores * (coneMosaicPatternSize*emPathLength*sizeOfDoubleInBytes + coneMosaicActivePatternSize*emPathLength*sizeOfDoubleInBytes/2)/(1024^3);
    trialBlockSize = floor(allowedRAMcompression*ramSizeGBytesAvailable/singleTrialMemoryGBytes);
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
            fprintf('trialBlock%d contains %d trials: %04d - %04d\n', iTrialBlock, numel(trialIndicesForThisBlock), firstTrial, lastTrial);
        end
        warndlg(...
            sprintf('CoresNum = %d; SystemRAM = %2.2fGB; estimated peak RAM = %2.2fGB', numberOfCores, ramSizeGBytes, trialBlockSize*singleTrialMemoryGBytes+ramUsedByOSGBytes), ...
            sprintf('%d trials computed in %d blocks, blockSize(1/last) = %d/%d', nTrials, numel(blockedTrialIndices), firstTrialBlockSize, lastTrialBlockSize));
    end
end
