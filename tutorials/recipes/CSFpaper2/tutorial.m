%% Tutorial listed in the Cottaris et al (2020) paper.
%% NPC, ISETBIO Team
%%
%% Generate a display for presenting stimuli and place it at a viewing distance of 57 cm
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 0.57);

%% Specify a Gabor stimulus 
stimParams = struct(...
    'spatialFrequencyCyclesPerDeg', 6, ...  % 6 cycles/deg
    'orientationDegs', 0, ...               % 0 degrees
    'widthDegs', 0.5, ...                   % 0.5 x 0.5 deg size
    'contrast', 0.8,...                     % 80% Michelson contrast
    'meanLuminanceCdPerM2', 40);            % 40 cd/m2 mean luminance

fprintf('\n1.Generating stimulus scene ...');
%% Generate an ISETBio scene describing this stimulus
stimulusScene = generateStimulusScene(stimParams, presentationDisplay);

%% Generate an ISETBio scene describing the background
stimParams.contrast = 0.0;
backgroundScene = generateStimulusScene(stimParams, presentationDisplay);


%% Realize the scenes into the particular LCD display
realizedStimulusScene = realizeSceneInDisplay(stimulusScene, presentationDisplay);
realizedBackgroundScene = realizeSceneInDisplay(backgroundScene, presentationDisplay);
fprintf(' Done !');


%% Generate wavefront-aberration derived human optics
fprintf('\n2.Generating optical image ...');
opticalImage = oiCreate('wvf human');

%% Compute the retinal images of the stimulus and the background scenes
stimulusOI = oiCompute(opticalImage, realizedStimulusScene);
backgroundOI = oiCompute(opticalImage, realizedBackgroundScene);
fprintf(' Done !');

%% Compute the stimulus temporal modulation function (square wave)
fprintf('\n3.Generating optical image sequence ...');
stimulusSamplingIntervalSeconds = 50/1000;     % 50 msec refresh time (20 Hz)
stimulusDurationSeconds = 150/1000;            % 150 msec duration
stimulusTimeAxisSeconds = -0.1:stimulusSamplingIntervalSeconds:0.2;   % Compute responses from -100ms to +200ms around the stimulus onset
stimONbins = stimulusTimeAxisSeconds>=0 & stimulusTimeAxisSeconds <= stimulusDurationSeconds-stimulusSamplingIntervalSeconds;
stimulusTemporalModulation = zeros(1, numel(stimulusTimeAxisSeconds));
stimulusTemporalModulation(stimONbins) = 1;

%% Compute the optical image sequence that simulates the stimulus presentation
theOIsequence = oiSequence(backgroundOI, stimulusOI, stimulusTimeAxisSeconds, stimulusTemporalModulation, 'composition', 'blend');

%% Visualize the sequence
%theOIsequence.visualize('montage'); 
fprintf(' Done !');

recomputeConeMosaic = false;
if (recomputeConeMosaic)
    %% Generate a hexagonal cone mosaic with ecc-based cone quantal efficiency
    fprintf('\n4.Generating cone mosaic ...');
    theConeMosaic = coneMosaicHex(7, ...             % hex lattice sampling factor
       'fovDegs', stimParams.widthDegs, ...       % match mosaic width to stimulus size
       'eccBasedConeDensity', true, ...           % cone density varies with eccentricity
       'eccBasedConeQuantalEfficiency', true, ... % cone quantal efficiency varies with eccentricity
       'integrationTime', 10/1000, ...            % 10 msec integration time
       'maxGridAdjustmentIterations', 50);        % terminate iterative lattice adjustment after 50 iterations
    save('coneMosaic.mat', 'theConeMosaic');
else
    fprintf('\n4.Loading cone mosaic ...');
    load('coneMosaic.mat');
    theConeMosaic = coneMosaic;
end
fprintf(' Done !');

%% Generate fixational eye movements for 1 trial and for the entire simulation time
fprintf('\n5.Generating fixational eye movements ...');
nTrials = 1;
eyeMovementsNum = ...
            theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
theEMPaths = theConeMosaic.emGenSequence(eyeMovementsNum, 'nTrials', nTrials);
fprintf(' Done !');

%% Response time axis
responseTimeAxis = (1:eyeMovementsNum)*theConeMosaic.integrationTime + theOIsequence.timeAxis(1);

%% Force eye movement to be at origin at t = 0
[~,idx] = min(abs(responseTimeAxis));
theEMPaths = bsxfun(@minus, theEMPaths, theEMPaths(:,idx,:));
theEMPathsDegs = theEMPaths * theConeMosaic.patternSampleSize(1) * 1e6 / theConeMosaic.micronsPerDegree;

%% Compute responses
fprintf('\n6.Computing responses ...');
for iTrial = 1:nTrials
        [absorptionsCountSequence(iTrial,:,:), photoCurrentSequence(iTrial,:,:)] = ...
                theConeMosaic.computeForOISequence(theOIsequence, ...
                'emPaths', theEMPaths(iTrial,:,:), ...
                'currentFlag', true);
end
fprintf(' Done !');

%% Visualize results
fprintf('\n7.Visualizing responses ...');
% Response visualization range
absorptionsRange = [0 max(absorptionsCountSequence(:))];
photoCurrentRange = prctile(photoCurrentSequence(:), [5 95]);
emRange = max(abs(theEMPathsDegs(:)))*[-1 1];
    
% Select some cones positioned around (0.05, -0.05)
horizontalPositionsDegs = 0.05+(-0.02:0.01:0.02);
verticalPositionsDegs = -0.05+(-0.02:0.01:0.02);
[xx,yy] = meshgrid(horizontalPositionsDegs, verticalPositionsDegs);
targetPosDegs = [xx(:) yy(:)];
[coneIndices, conePositionsDegs, coneTypes, coneIndicesInSerializedList] = ...
            theConeMosaic.indicesForConesAtPositions(targetPosDegs);
        
lConeIndices = coneIndicesInSerializedList(coneTypes == 2);
mConeIndices = coneIndicesInSerializedList(coneTypes == 3);
sConeIndices = coneIndicesInSerializedList(coneTypes == 4);

        
        
visualizedTrial = 1;
visualizedAbsorptions = squeeze(absorptionsCountSequence(visualizedTrial,:,:));
visualizedPhotocurrents = squeeze(photoCurrentSequence(visualizedTrial,:,:));
visualizedEMPath = squeeze(theEMPathsDegs(visualizedTrial,:,:));

% Compute mean responses during the pre-stimulus interval
nonCausalTimes = find(responseTimeAxis<0);
meanAbsorptions = mean(visualizedAbsorptions(:, nonCausalTimes), 2);
meanPhotocurrents = mean(visualizedPhotocurrents(:, nonCausalTimes), 2);
% Subtract mean responses to compute differential responses (modulations)
differentialAbsorptions = bsxfun(@minus, visualizedAbsorptions, meanAbsorptions);
differentialPhotocurrents = bsxfun(@minus, visualizedPhotocurrents, meanPhotocurrents);
diffAbsorptionsRange = max(abs(differentialAbsorptions(:)))*[-1 1];
diffPhotoCurrentRange = max(abs(differentialPhotocurrents(:)))*[-1 1];

% Generate movie
hFig = figure(1); clf;
set(hFig, 'Position', [10 10 900 900]);
videoOBJ = VideoWriter('responseVideo', 'MPEG-4'); % H264 format
videoOBJ.FrameRate = 30;
videoOBJ.Quality = 100;
videoOBJ.open();
        
% Plot cone mosaic
axHandle = subplot(3,2,1);
theConeMosaic.visualizeGrid('axesHandle', axHandle, ...
        'apertureShape', 'hexagons', 'ticksInVisualDegs', true, 'tickInc', 0.1);
    

for tBin = 2:numel(responseTimeAxis)
    % Plot fixational eye movement trajectory
    fprintf('\n\tFrame %d of %d ...', tBin, numel(responseTimeAxis));
    subplot(3,2,2);
    plot(visualizedEMPath(1:tBin,1), -visualizedEMPath(1:tBin,2), 'b-', 'LineWidth', 1.5); hold on;
    plot(visualizedEMPath(tBin,1), -visualizedEMPath(tBin,2), 'rs', 'LineWidth', 1.5); hold off;
    set(gca, 'XLim', [-0.2 0.2], 'YLim', [-0.2 0.2], 'XTick', -0.2:0.1:0.2, 'YTick', -0.2:0.1:0.2,  'FontSize', 18);
    grid on; box on;
    xlabel('space (degs)');
    ylabel('space (degs)');
    title(sprintf('%2.0f ms', responseTimeAxis(tBin)*1000));
    axis 'square';

    % Plot cone mosaic excitation map (absorptions)
    axHandle = subplot(3,2,3);
    theConeMosaic.renderActivationMap(axHandle, differentialAbsorptions(:,tBin), ...
        'mapType', 'modulated hexagons', 'signalRange', diffAbsorptionsRange, ...
        'tickInc', 0.1);

    % Plot cone mosaic excitation map (photocurrents)
    axHandle = subplot(3,2,5);
    theConeMosaic.renderActivationMap(axHandle, differentialPhotocurrents(:,tBin), ...
        'mapType', 'modulated hexagons', 'signalRange', diffPhotoCurrentRange, ...
        'tickInc', 0.1);

    % Plot time series of absorptions
    subplot(3,2,4);
    if (~isempty(lConeIndices))
        plot(responseTimeAxis(1:tBin)*1000, differentialAbsorptions(lConeIndices, 1:tBin), 'r-', 'LineWidth', 1.5); 
    end
    hold on;
    if (~isempty(mConeIndices))
        plot(responseTimeAxis(1:tBin)*1000, differentialAbsorptions(mConeIndices, 1:tBin), 'g-', 'LineWidth', 1.5); 
    end
    if (~isempty(sConeIndices))
        plot(responseTimeAxis(1:tBin)*1000, differentialAbsorptions(sConeIndices, 1:tBin), 'b-', 'LineWidth', 1.5); 
    end
    hold off

    set(gca, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)]*1000, 'YLim', diffAbsorptionsRange);
    set(gca, 'FontSize', 16);
    xlabel('time (ms)');
    ylabel('isomerization modulation (R*)');

    % Plot time series of photocurrents
    subplot(3,2,6);
    if (~isempty(lConeIndices))
        plot(responseTimeAxis(1:tBin)*1000, differentialPhotocurrents(lConeIndices, 1:tBin), 'r-', 'LineWidth', 1.5); 
    end
    hold on
    if (~isempty(mConeIndices))
        plot(responseTimeAxis(1:tBin)*1000, differentialPhotocurrents(mConeIndices, 1:tBin), 'g-','LineWidth', 1.5);
    end
    if (~isempty(sConeIndices))
        plot(responseTimeAxis(1:tBin)*1000, differentialPhotocurrents(sConeIndices, 1:tBin), 'b-','LineWidth', 1.5);
    end
    hold off
    set(gca, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)]*1000, 'YLim', diffPhotoCurrentRange );
    set(gca, 'FontSize', 16);
    xlabel('time (ms)');
    ylabel('pCurrent modulation (pAmps)');

    drawnow;
    % Add video frame
    videoOBJ.writeVideo(getframe(hFig));
end

% Close video stream
videoOBJ.close();

    