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
%   'trialBlocks' - How many blocks to split the testDirectionParams.trialsNum into. Default: 1 (no blocking). 
%               This only has an effect with @coneMosaicHex mosaics and when nTrials>1 and it is useful with 
%               large mosaics x lots of trials, in which case the absorptions matrix does not fit in the RAM.

%  'useSinglePrecision' - true/false (default true) use single precision to represent isomerizations and photocurrent

% 7/10/16  npc Wrote it.

%% Parse arguments
p = inputParser;
p.addParameter('seed',1, @isnumeric);   
p.addParameter('workerID', [], @isnumeric);
p.addParameter('useSinglePrecision',true,@islogical);
p.addParameter('trialBlocks', 1, @isnumeric);
p.parse(varargin{:});
currentSeed = p.Results.seed;
trialBlocks = p.Results.trialBlocks;

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
fudgeFactor = 5;
[numberOfCores, ramSizeGBytes, sizeOfDoubleInBytes] = determineSystemResources();
maxMemoryRequiredGBytes = fudgeFactor * 2*(numberOfCores * numel(theMosaic.pattern) * nTrials * sizeOfDoubleInBytes)/(1024^3);
desiredTrialsPerBlock = floor(min([nTrials nTrials / (maxMemoryRequiredGBytes/ramSizeGBytes)]));
trialBlocksForParforLoop = ceil(nTrials / desiredTrialsPerBlock);
%warndlg(...
%    sprintf('CoresNum = %d; SystemRAM = %2.2fGB; requiredRAM (all trials) = %2.2fGB', numberOfCores, ramSizeGBytes, %maxMemoryRequiredGBytes), ...
%    sprintf('nTrials = %d, trialBlocksForParforLoop = %d', nTrials, trialBlocksForParforLoop));
    
if (trialBlocks == -1)
    trialBlocks = trialBlocksForParforLoop;
end

[isomerizations, photocurrents] = ...
     theMosaic.computeForOISequence(theOIsequence, ...
                    'seed', currentSeed, ...
                    'emPaths', theEMpaths, ...
                    'trialBlocks', trialBlocks, ...
                    'currentFlag', true, ...
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
