function hFig = visualizeResponseInstance(conditionDir, responseInstance, responseLabel, theMosaic, iTrial, trialsNum, renderVideo)
% hFig = visualizeResponseInstance(conditionDir, responseInstance, theMosaic, iTrial, trialsNum)
% 
% Visualize the central (in time) part of a response instance
%
%  7/9/16  npc Wrote it.
%
    
    % Retrieve isomerization rate, photocurrent, and eye movement sequence
    isomerizationRate = responseInstance.theMosaicIsomerizations/theMosaic.integrationTime;
    photocurrentSequence = responseInstance.theMosaicPhotoCurrents;
    eyeMovementSequence = responseInstance.theMosaicEyeMovements;
    timeAxis = responseInstance.timeAxis;
    
    % Determine time bins to visualize. Only up to max response here.
    tPos = find(timeAxis > 0);
    timeStepsVisualized = 1:tPos(1);
    
    % Determine plotting ranges
    isomerizationRateRange =  [min(isomerizationRate(:)) max(isomerizationRate(:))]; % prctile(isomerizationRate(:), [1 99]); % prctile(isomerizationRate(:),99) + [-6000 0];    % [4500 10500]; % prctile(isomerizationRate(:), [1 99])
    photocurrentRange =  [min(photocurrentSequence(:)) max(photocurrentSequence(:))]; % prctile(photocurrentSequence(:), [1 99]); % prctile(photocurrentSequence(:), 1) + [0 20];      % [-75 -55];    % prctile(photocurrent(:), [1 99]);
    
    % Compute mosaic spatial axes
    coneRows = size(isomerizationRate,1);
    coneCols = size(isomerizationRate,2);
    mosaicXaxis = linspace(-coneCols/2, coneCols/2, coneCols);
    mosaicYaxis = linspace(-coneRows/2, coneRows/2, coneRows);

    
    hFig = figure(10); 
    set(hFig, 'Name', sprintf('%s [Trial %d / %d]', responseLabel, iTrial, trialsNum), 'Position', [10 10 1070 520], 'Color', [1 1 1]);
    clf; colormap(gray(1024));

    if (renderVideo)
        % Open video stream
        videoDir = colorGaborDetectOutputDir(conditionDir,'videos');
        if (~exist(fullfile(videoDir,responseLabel),'dir'))
            mkdir(fullfile(videoDir,responseLabel));
        end
        videoFilename = fullfile(videoDir,responseLabel,sprintf('Instance%dMovie.m4v', iTrial));
        writerObj = VideoWriter(videoFilename, 'MPEG-4'); % H264 format
        writerObj.FrameRate = 15; 
        writerObj.Quality = 50;
        writerObj.open();
    end
    
    for timeStep = timeStepsVisualized
        subplot('Position', [0.01 0.03 0.43 0.94]);
        imagesc(mosaicXaxis, mosaicYaxis, isomerizationRate(:,:,timeStep), isomerizationRateRange);
        hold on;
        idx = [max([1 timeStep-10]) : timeStep];
        plot(eyeMovementSequence(idx,1), -eyeMovementSequence(idx,2), 'w-', 'Color', [1.0 0.5 0.5], 'LineWidth', 4.0);
        plot(eyeMovementSequence(idx,1), -eyeMovementSequence(idx,2), 'r.-', 'LineWidth', 2.0);
        hold off;
        axis 'image'; axis 'xy'
        %xlabel(sprintf('%2.0f microns (%2.2f deg)', theMosaic.width*1e6, theMosaic.fov(1)), 'FontSize', 14, 'FontName', 'Menlo');
        set(gca, 'CLim', isomerizationRateRange, 'XTick', [], 'YTick', []);
        hCbar = colorbar(); % 'Ticks', cbarStruct.ticks, 'TickLabels', cbarStruct.tickLabels);
        hCbar.Orientation = 'vertical'; 
        hCbar.Label.String = 'isomerization rate (R*/cone/sec)'; 
        hCbar.FontSize = 14; 
        hCbar.FontName = 'Menlo'; 
        hCbar.Color = [0.2 0.2 0.2];
        title(sprintf('%s\nisomerization rate map (t=%2.2f ms)', responseLabel, timeAxis(timeStep)*1000), 'FontSize', 16, 'FontName', 'Menlo');

        subplot('Position', [0.53 0.03 0.43 0.94]);
        imagesc(photocurrentSequence(:,:,timeStep));
        % xlabel(sprintf('%2.0f microns (%2.2f deg)', theMosaic.width*1e6, theMosaic.fov(1)), 'FontSize', 14, 'FontName', 'Menlo');
        axis 'image'; axis 'xy'
        set(gca, 'CLim', photocurrentRange, 'XTick', [], 'YTick', []);
        hCbar = colorbar(); % 'Ticks', cbarStruct.ticks, 'TickLabels', cbarStruct.tickLabels);
        hCbar.Orientation = 'vertical'; 
        hCbar.Label.String = 'photocurrent (pAmps)'; 
        hCbar.FontSize = 14; 
        hCbar.FontName = 'Menlo'; 
        hCbar.Color = [0.2 0.2 0.2];
        title(sprintf('%s\nphotocurrent map (t=%2.2f ms)', responseLabel, timeAxis(timeStep)*1000), 'FontSize', 16, 'FontName', 'Menlo');
        drawnow;
        if (renderVideo)
            writerObj.writeVideo(getframe(hFig));
        end
    end
    
    if (renderVideo)
        writerObj.close();
        fprintf('Movie saved in %s\n', videoFilename);
    end
    
end
