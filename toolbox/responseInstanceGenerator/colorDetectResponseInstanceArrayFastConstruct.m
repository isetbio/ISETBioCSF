function [responseStruct, osImpulseResponseFunctions, osImpulseReponseFunctionTimeAxis, osMeanCurrents] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic, varargin)
% [responseStruct, osImpulseResponseFunctions, osMeanCurrents] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, nTrials, spatialParams, backgroundParams, colorModulationParams, temporalParams, theOI, theMosaic)
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
%  'displayTrialBlockPartitionDiagnostics' -true/false (default false). Wether to report time elapsed for each block
%  'centeredEMPaths' - true/false (default false) 
%               Controls wether the eye movement paths start at (0,0) (default) or wether they are centered around (0,0)
%  'osImpulseResponseFunctions' - the LMS impulse response filters to be used (default: [])
%  'osMeanCurrents' - the steady-state LMS currents caused by the mean absorption LMS rates
%  'computePhotocurrentResponseInstances' - true/false (default true) wether to compute photocurrent response instances
%  'computeNoiseFreeSignals' - true/false (default true) wether to compute the noise free isomerizations and photocurrents
%  'useSinglePrecision' - true/false (default true) use single precision to represent isomerizations and photocurrent
%  'visualizeSpatialScheme' - true/false (default false) wether to co-visualize the optical image and the cone mosaic
%  'paramsList' - the paramsList associated for this direction (used for exporting response figures to the right directory)
% 7/10/16  npc Wrote it.

%% Parse arguments
p = inputParser;
p.addParameter('seed',1, @isnumeric);   
p.addParameter('workerID', [], @isnumeric);
p.addParameter('displayTrialBlockPartitionDiagnostics', false, @islogical);
p.addParameter('useSinglePrecision',true,@islogical);
p.addParameter('centeredEMPaths',false, @islogical); 
p.addParameter('osImpulseResponseFunctions', [], @isnumeric);
p.addParameter('osMeanCurrents', [], @isnumeric);
p.addParameter('computePhotocurrentResponseInstances', true, @islogical);
p.addParameter('computeNoiseFreeSignals', true, @islogical);
p.addParameter('visualizeSpatialScheme', false, @islogical);
p.addParameter('visualizeOIsequence', false, @islogical);
p.addParameter('visualizeMosaicWithFirstEMpath', false, @islogical);
p.addParameter('paramsList', {}, @iscell);
p.parse(varargin{:});
currentSeed = p.Results.seed;


%% Start computation time measurement
if (p.Results.displayTrialBlockPartitionDiagnostics)
    tic
end

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
if (isnan(temporalParams.windowTauInSeconds))
    [stimulusTimeAxis, stimulusModulationFunction, ~] = squareTemporalWindowCreate(temporalParams);
else
    [stimulusTimeAxis, stimulusModulationFunction, ~] = gaussianTemporalWindowCreate(temporalParams);
end

%% Compute the oiSequence
if (strcmp(spatialParams.spatialType, 'pedestalDisk'))
    theOIsequence = oiSequence(oiBackground, oiModulated, stimulusTimeAxis, stimulusModulationFunction, ...
        'composition', 'xor');
else
    theOIsequence = oiSequence(oiBackground, oiModulated, stimulusTimeAxis, stimulusModulationFunction, ...
        'composition', 'blend');
end


%%  Visualize the oiSequence
if (p.Results.visualizeOIsequence)
    [~, hFig] = theOIsequence.visualize('format', 'montage', 'showIlluminanceMap', true);
    % Save figure
    theProgram = mfilename;
    rwObject = IBIOColorDetectReadWriteBasic;
    data = 0;
    fileName = sprintf('OISequence');
    rwObject.write(fileName, data, p.Results.paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
       
    % visualize the scene and the optical image
    visualizedPercentange = 0.2;
    sceneLuminance = sceneGet(modulatedScene, 'luminance');
    spatialSupport = sceneGet(modulatedScene, 'spatial support');
    xSupport = squeeze(spatialSupport(1,:,1));
    ySupport = squeeze(spatialSupport(:,1,2));
    idx = find((abs(squeeze(spatialSupport(:,:,1))) < max(xSupport)*visualizedPercentange) & ...
               (abs(squeeze(spatialSupport(:,:,2))) < max(ySupport)*visualizedPercentange));
    maxSceneLuminance = max(max(sceneLuminance(idx)));
    minSceneLuminance = min(min(sceneLuminance(idx)));      
    
    oiIlluminance = oiGet(oiModulated,'illuminance');
	oiSpatialSupport = oiGet(oiModulated, 'spatial support');
    oiXSupport = squeeze(oiSpatialSupport(1,:,1));
    oiYSupport = squeeze(oiSpatialSupport(:,1,2));
    idx = find((abs(squeeze(oiSpatialSupport(:,:,1))) < max(oiXSupport)*visualizedPercentange) & ...
               (abs(squeeze(oiSpatialSupport(:,:,2))) < max(oiYSupport)*visualizedPercentange));
    maxOIILuminance = max(max(oiIlluminance(idx)));
    minOIILuminance = min(min(oiIlluminance(idx)));
    
    
    hFig = figure(99); clf;
    set(hFig, 'Position', [10 10 1000 700], 'Color', [1 1 1]);
    subplot(1,2,1)
    pcolor(xSupport, ySupport, (sceneLuminance-minSceneLuminance)/(maxSceneLuminance-minSceneLuminance));
    set(gca, 'CLim', [0 1])
    hold on;
    plot([0 0], [ySupport(1) ySupport(end)], 'w-');
    plot([xSupport(1) xSupport(end)], [0 0], 'w-');
    hold off
    axis 'image';
    set(gca, 'XLim', [-max(xSupport) max(xSupport)]*visualizedPercentange, 'XTick', [0], ...
             'YLim', [-max(ySupport) max(ySupport)]*visualizedPercentange, 'YTick', [0]);
    set(gca, 'FontSize', 16);
    title('scene lluminance');
    subplot(1,2,2)
    pcolor(oiXSupport, oiYSupport, (oiIlluminance-minOIILuminance)/(maxOIILuminance-minOIILuminance));
    set(gca, 'CLim', [0 1])
    hold on;
    plot([0 0], [ySupport(1) ySupport(end)], 'w-');
    plot([xSupport(1) xSupport(end)], [0 0], 'w-');
    hold off
    axis 'image';
    set(gca, 'XLim', [-max(oiXSupport) max(oiXSupport)]*visualizedPercentange, 'XTick', [0], ...
             'YLim', [-max(oiYSupport) max(oiYSupport)]*visualizedPercentange, 'YTick', [0]);
    set(gca, 'FontSize', 16);
    title('optical image illuminance');
    colormap(jet(1024));
    drawnow
    
    fileName = sprintf('SceneLuminanceAndIlluminance');
    rwObject.write(fileName, data, p.Results.paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
end

%% Co-visualize the optical image and the cone mosaic
if (p.Results.visualizeSpatialScheme)
    [~, timeBinOfPeakModulation] = max(abs(stimulusModulationFunction));
    thePeakOI = theOIsequence.frameAtIndex(timeBinOfPeakModulation);
    visualizeStimulusAndConeMosaic(theMosaic, thePeakOI, p.Results.paramsList);
end
    
% Clear oiModulated and oiBackground to save space
varsToClear = {'oiBackground', 'oiModulated'};
clear(varsToClear{:});

%% Generate eye movement paths for all instances
eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theMosaic.integrationTime);
[theEMpaths, theEMpathsMicrons] = colorDetectMultiTrialEMPathGenerate(...
    theMosaic, nTrials, eyeMovementsNum, temporalParams.emPathType, ...
    'centeredEMPaths', p.Results.centeredEMPaths, ...
    'seed', currentSeed);

if (p.Results.visualizeMosaicWithFirstEMpath)
    visualizedTrial = 1;
    hFig = theMosaic.visualizeGrid(...
        'overlayEMpathMicrons', squeeze(theEMpathsMicrons(visualizedTrial,:,:)), ...  % expects emPath in microns, not cone units
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'lightCollectingArea', ...
        'labelConeTypes', true,...
        'generateNewFigure', true);
    % Save figure
    theProgram = mfilename;
    rwObject = IBIOColorDetectReadWriteBasic;
    data = 0;
    fileName = sprintf('MosaicAndFirstEMpath');
    rwObject.write(fileName, data, p.Results.paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
end

[isomerizations, photocurrents, osImpulseResponseFunctions, osMeanCurrents] = ...
     theMosaic.computeForOISequence(theOIsequence, ...
                    'seed', currentSeed, ...
                    'emPaths', theEMpaths, ...
                    'interpFilters', p.Results.osImpulseResponseFunctions, ...
                    'meanCur', p.Results.osMeanCurrents, ...
                    'trialBlockSize', nTrials, ...      % Do all trials in a single block
                    'currentFlag', p.Results.computePhotocurrentResponseInstances, ...
                    'workDescription', sprintf('%d trials of noisy responses', nTrials),...
                    'workerID', p.Results.workerID);
                
if (isempty(p.Results.osImpulseResponseFunctions))
    % Since we compute the impulse response functions, the 
    % mosaic.interpFilterTimeAxis has already being computed for us
    osImpulseReponseFunctionTimeAxis = theMosaic.interpFilterTimeAxis;
else
    % Since we do not compute the impulse response functions, the
    % mosaic.interpFilterTimeAxis is not computed, so compute it here
    osImpulseReponseFunctionTimeAxis = (0:(size(p.Results.osImpulseResponseFunctions,1)-1))*theMosaic.integrationTime;
end

isomerizationsTimeAxis = theMosaic.timeAxis + theOIsequence.timeAxis(1);
photoCurrentTimeAxis = isomerizationsTimeAxis;
 
%% Remove unwanted portions of the responses
isomerizationsTimeIndicesToKeep = find(abs(isomerizationsTimeAxis-temporalParams.secondsToIncludeOffset) <= temporalParams.secondsToInclude/2);
photocurrentsTimeIndicesToKeep = find(abs(photoCurrentTimeAxis-temporalParams.secondsToIncludeOffset) <= temporalParams.secondsToInclude/2);
if (isa(theMosaic , 'coneMosaicHex'))
     isomerizations = isomerizations(:,:,isomerizationsTimeIndicesToKeep);
     if (~isempty(photocurrents))
        photocurrents = photocurrents(:,:,photocurrentsTimeIndicesToKeep);
     end
else
     isomerizations = isomerizations(:,:,:,isomerizationsTimeIndicesToKeep);
     if (~isempty(photocurrents))
        photocurrents = photocurrents(:,:,:,photocurrentsTimeIndicesToKeep);
     end
end

%% Form the responseInstanceArray struct
if (p.Results.useSinglePrecision)
    responseStruct.responseInstanceArray.theMosaicIsomerizations = single(isomerizations);
    responseStruct.responseInstanceArray.theMosaicPhotocurrents = single(photocurrents);
else
    responseStruct.responseInstanceArray.theMosaicIsomerizations = isomerizations;
    responseStruct.responseInstanceArray.theMosaicPhotocurrents = photocurrents;
end
responseStruct.responseInstanceArray.theMosaicEyeMovements = theEMpaths(:,isomerizationsTimeIndicesToKeep,:);
responseStruct.responseInstanceArray.timeAxis = isomerizationsTimeAxis(isomerizationsTimeIndicesToKeep);
responseStruct.responseInstanceArray.photocurrentTimeAxis = photoCurrentTimeAxis(photocurrentsTimeIndicesToKeep);

if (~p.Results.computeNoiseFreeSignals)
    responseStruct.noiseFreeIsomerizations = [];
    responseStruct.noiseFreePhotocurrents = [];
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
originalPhotocurrentNoiseFlag = theMosaic.os.noiseFlag;

% Set noiseFlags to none
theMosaic.noiseFlag = 'none';
theMosaic.os.noiseFlag = 'none';

% Compute the noise-free responses
[responseStruct.noiseFreeIsomerizations, responseStruct.noiseFreePhotocurrents] = ...
     theMosaic.computeForOISequence(theOIsequence, ...
                    'emPaths',theEMpaths, ...
                    'interpFilters', p.Results.osImpulseResponseFunctions, ...
                    'meanCur', p.Results.osMeanCurrents, ...
                    'currentFlag', true, ...
                    'workDescription', 'noise-free responses', ...
                    'workerID', p.Results.workerID);

%% Remove unwanted portions of the noise-free responses
if (isa(theMosaic , 'coneMosaicHex'))
     responseStruct.noiseFreeIsomerizations = responseStruct.noiseFreeIsomerizations(:,:,isomerizationsTimeIndicesToKeep);
     responseStruct.noiseFreePhotocurrents = responseStruct.noiseFreePhotocurrents(:,:,photocurrentsTimeIndicesToKeep);
else
     responseStruct.noiseFreeIsomerizations = responseStruct.noiseFreeIsomerizations(:,:,:,isomerizationsTimeIndicesToKeep);
     responseStruct.noiseFreePhotocurrents = responseStruct.noiseFreePhotocurrents(:,:,:,photocurrentsTimeIndicesToKeep);
end

% Restore noiseFlags to none
theMosaic.noiseFlag = originalIsomerizationNoiseFlag;
theMosaic.os.noiseFlag = originalPhotocurrentNoiseFlag;

% Store
responseStruct.noiseFreeIsomerizations = squeeze(responseStruct.noiseFreeIsomerizations);
if (~isempty(responseStruct.noiseFreePhotocurrents))
    responseStruct.noiseFreePhotocurrents = squeeze(responseStruct.noiseFreePhotocurrents);
end

%% Report time taken
if (p.Results.displayTrialBlockPartitionDiagnostics)
    if isempty(p.Results.workerID)
        fprintf('<strong> Response instance array computation (%d instances) took %2.3f minutes. </strong> \n', nTrials, toc/60);
    else
        fprintf('<strong> Response instance array computation (%d instances) in worker %d took %2.3f minutes. </strong> \n', nTrials, p.Results.workerID, toc/60);
    end
end

end