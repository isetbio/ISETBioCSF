%% Tutorial listed in the Cottaris et al (2019B) paper.
%% NPC, ISETBIO Team
%%
%% Generate a display for presenting stimuli and place it at a viewing distance of 57 cm
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 0.57);

%% Specify a Gabor stimulus 
stimParams = struct(...
    'spatialFrequencyCyclesPerDeg', 4, ...  % 4 cycles/deg
    'orientationDegs', 0, ...               % 0 degrees
    'widthDegs', 0.4, ...                   % 0.x4 x 0.4 size
    'contrast', 0.9,...                     % 80% Michelson contrast
    'meanLuminanceCdPerM2', 40);            % 40 cd/m2 mean luminance

%% Generate an ISETBio scene describing this stimulus
stimulusScene = generateStimulusScene(stimParams, presentationDisplay);

%% Generate an ISETBio scene describing the background
stimParams.contrast = 0.0;
backgroundScene = generateStimulusScene(stimParams, presentationDisplay);


%% Realize the scenes into the particular LCD display
realizedStimulusScene = realizeSceneInDisplay(stimulusScene, presentationDisplay);
realizedBackgroundScene = realizeSceneInDisplay(backgroundScene, presentationDisplay);

%% Generate wavefront-aberration derived human optics
opticalImage = oiCreate('wvf human');

%% Compute the retinal images of the stimulus and the background scenes
stimulusOI = oiCompute(opticalImage, realizedStimulusScene);
backgroundOI = oiCompute(opticalImage, realizedBackgroundScene);

%% Compute the stimulus temporal modulation function (square wave)
stimulusSamplingIntervalSeconds = 50/1000;     % 50 msec refresh time (20 Hz)
stimulusDurationSeconds = 150/1000;            % 150 msec duration
stimulusTimeAxisSeconds = -0.1:stimulusSamplingIntervalSeconds:0.2;   % Compute responses from -100ms to +200ms around the stimulus onset
stimONbins = stimulusTimeAxisSeconds>=0 & stimulusTimeAxisSeconds <= stimulusDurationSeconds-stimulusSamplingIntervalSeconds;
stimulusTemporalModulation = zeros(1, numel(stimulusTimeAxisSeconds));
stimulusTemporalModulation(stimONbins) = 1;
    
%% Compute the optical image sequence that simulates the stimulus presentation
theOIsequence = oiSequence(backgroundOI, stimulusOI, stimulusTimeAxisSeconds, stimulusTemporalModulation, 'composition', 'blend');
theOIsequence.visualize('montage'); 


recomputeConeMosaic = false;
if (recomputeConeMosaic)
    %% Generate a hexagonal cone mosaic with ecc-based cone quantal efficiency
    theConeMosaic = coneMosaicHex(7, ...             % hex lattice sampling factor
       'fovDegs', stimParams.widthDegs, ...       % match mosaic width to stimulus size
       'eccBasedConeDensity', true, ...           % cone density varies with eccentricity
       'eccBasedConeQuantalEfficiency', true, ... % cone quantal efficiency varies with eccentricity
       'integrationTime', 10/1000, ...            % 10 msec integration time
       'maxGridAdjustmentIterations', 50);        % terminate iterative lattice adjustment after 50 iterations
    save('coneMosaic.mat', 'theConeMosaic');
else
    load('coneMosaic.mat');
    theConeMosaic = coneMosaic;
end

%% Generate fixational eye movements for 1 trial and for the entire simulation time
nTrials = 1;
eyeMovementsNum = ...
            theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
theEMPaths = theConeMosaic.emGenSequence(eyeMovementsNum, 'nTrials', nTrials);


%% Response time axis
responseTimeAxis = (1:eyeMovementsNum)*theConeMosaic.integrationTime + theOIsequence.timeAxis(1);

%% Force eye movement to be at origin at t = 0
[~,idx] = min(abs(responseTimeAxis));
theEMPaths = bsxfun(@minus, theEMPaths, theEMPaths(:,idx,:));
    
%% Compute responses
for iTrial = 1:nTrials
        [absorptionsCountSequence(iTrial,:,:), photoCurrentSequence(iTrial,:,:)] = ...
                theConeMosaic.computeForOISequence(theOIsequence, ...
                'emPaths', theEMPaths(iTrial,:,:), ...
                'currentFlag', true);
end
    
%% Visualize results

% Response visualization range
absorptionsRange = [0 max(absorptionsCountSequence(:))];
photoCurrentRange = prctile(photoCurrentSequence(:), [5 95]);
emRange = max(abs(theEMPaths(:)))*[-1 1];
    
% Indices of L/M/S cones
    lConeIndices = find(theConeMosaic.coneTypesHexGrid == 2);
    mConeIndices = find(theConeMosaic.coneTypesHexGrid == 3);
    sConeIndices = find(theConeMosaic.coneTypesHexGrid == 4);
    
visualizedTrial = 1;
visualizedAbsorptions = squeeze(absorptionsCountSequence(visualizedTrial,:,:));
visualizedPhotocurrents = squeeze(photoCurrentSequence(visualizedTrial,:,:));
visualizedEMPath = squeeze(theEMPaths(visualizedTrial,:,:));
    
hFig = figure(1); clf;
videoOBJ = VideoWriter('responseVideo', 'MPEG-4'); % H264 format
videoOBJ.FrameRate = 30;
videoOBJ.Quality = 100;
videoOBJ.open();
        
for tBin = 1:numel(responseTimeAxis)

    subplot(3,2,1);
    plot(visualizedEMPath(1:tBin,1), -visualizedEMPath(1:tBin,2), 'k-'); hold on;
    plot(visualizedEMPath(tBin,1), -visualizedEMPath(tBin,2), 'rs'); hold off;
    set(gca, 'XLim', emRange, 'YLim', emRange);
    title(sprintf('%2.0f ms', responseTimeAxis(tBin)*1000));
    axis 'square';

    axHandle = subplot(3,2,3);
    theConeMosaic.renderActivationMap(axHandle, visualizedAbsorptions(:,tBin), ...
        'mapType', 'modulated disks', 'signalRange', absorptionsRange);

    axHandle = subplot(3,2,5);
    theConeMosaic.renderActivationMap(axHandle, visualizedPhotocurrents(:,tBin), ...
        'mapType', 'modulated disks', 'signalRange', photoCurrentRange);

    subplot(3,2,4);

    plot(responseTimeAxis(1:tBin)*1000, visualizedAbsorptions(lConeIndices, 1:tBin), 'r-'); hold on
    plot(responseTimeAxis(1:tBin)*1000, visualizedAbsorptions(mConeIndices, 1:tBin), 'g-');
    plot(responseTimeAxis(1:tBin)*1000, visualizedAbsorptions(sConeIndices, 1:tBin), 'b-');
    hold off

    set(gca, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)]*1000, 'YLim', [min(visualizedAbsorptions(:)) max(visualizedAbsorptions(:))]);
    set(gca, 'FontSize', 16);
    xlabel('time (ms)');
    ylabel('isomerizations');

    subplot(3,2,6);
    plot(responseTimeAxis(1:tBin)*1000, visualizedPhotocurrents(lConeIndices, 1:tBin), 'r-'); hold on
    plot(responseTimeAxis(1:tBin)*1000, visualizedPhotocurrents(mConeIndices, 1:tBin), 'g-');
    plot(responseTimeAxis(1:tBin)*1000, visualizedPhotocurrents(sConeIndices, 1:tBin), 'b-');
    hold off
    set(gca, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)]*1000, 'YLim', [min(visualizedPhotocurrents(:)) max(visualizedPhotocurrents(:))]);
    set(gca, 'FontSize', 16);
    xlabel('time (ms)');
    ylabel('pCurrent');

    drawnow;
    % Add video frame
    videoOBJ.writeVideo(getframe(hFig));
end

% Close video stream
videoOBJ.close();
    