function visualizeMosaicResponseSequence(conditionDir, signalName, mosaicResponseSequence, eyeMovementSequence, coneTypes, timeAxis, mosaicSize, mosaicFOV, integrationTimeInSeconds, movieFileName)
% visualizeMosaicResponseSequence(conditionDir, signalName, mosaicResponseSequence, eyeMovementSequence, coneTypes, timeAxis, mosaicSize, mosaicFOV, movieFileName)
%
% Visualize the time cource of a mosaic response (and possibly the eye movement path)
% and generate a video of it.
%
%
%  7/12/16  npc Wrote it.
%
    
    % Determine ranges for plotting
    mosaicRows = size(mosaicResponseSequence,1);
    mosaicCols = size(mosaicResponseSequence,2);
    mosaicXaxis = linspace(-mosaicCols/2, mosaicCols/2, mosaicCols);
    mosaicYaxis = linspace(-mosaicRows/2, mosaicRows/2, mosaicRows);
    
    % Open video stream
    videoDir = colorGaborDetectOutputDir(conditionDir,'videos');
    videoFilename = fullfile(videoDir, sprintf('%s.m4v', movieFileName));
    writerObj = VideoWriter(videoFilename, 'MPEG-4'); % H264 format
    writerObj.FrameRate = 15; 
    writerObj.Quality = 100;
    writerObj.open();
    
    
    hFig = figure(11); clf;
    set(hFig, 'Position', [10 10 1240 445], 'Color', [1 1 1]);
    
    
    % Make cone mask to separate submosaic responses
    coneMask = zeros(mosaicRows, mosaicCols, 3);
    for coneType = 1 : 3 % L, M, S
        coneMask(:,:,coneType) = double(coneTypes == coneType+1);
    end
    
    responseRange = prctile(mosaicResponseSequence(:), [1 99]);
    colorbarTicks = linspace(1,1023, 10);
    colorbarTickLabels = linspace(responseRange(1), responseRange(2), 10);
    
    normalizedMosaicResponseSequence = (mosaicResponseSequence-responseRange(1))/(responseRange(2) - responseRange(1));
    normalizedMosaicResponseSequence(normalizedMosaicResponseSequence<0) = 0;
    normalizedMosaicResponseSequence(normalizedMosaicResponseSequence>1) = 1;
    normalizedMosaicResponseSequence = 1+round(1022*normalizedMosaicResponseSequence);
    
    colorLUT = (jet(1024)).^0.8;
    colorLUT(1,:) = [0 0 0];
    colormap(colorLUT);
    
    LconeMosaicResponseSequence = normalizedMosaicResponseSequence .* repmat(squeeze(coneMask(:,:,1)), [1 1 size(mosaicResponseSequence,3)]);
    MconeMosaicResponseSequence = normalizedMosaicResponseSequence .* repmat(squeeze(coneMask(:,:,2)), [1 1 size(mosaicResponseSequence,3)]);
    SconeMosaicResponseSequence = normalizedMosaicResponseSequence .* repmat(squeeze(coneMask(:,:,3)), [1 1 size(mosaicResponseSequence,3)]);

    
    fprintf('Rendering responses ...\n');
    for timeStep = 1:size(mosaicResponseSequence,3)
        
        LconeMosaicResponse = squeeze(LconeMosaicResponseSequence(:,:,timeStep));
        MconeMosaicResponse = squeeze(MconeMosaicResponseSequence(:,:,timeStep));
        SconeMosaicResponse = squeeze(SconeMosaicResponseSequence(:,:,timeStep));
        
        % L-cone mosaic subplot
        subplot('Position', [0.04 0.09 0.28 0.82]);
        imagesc(mosaicXaxis, mosaicYaxis, LconeMosaicResponse);
        % Superimpose the eye movement path to this point
        if (~isempty(eyeMovementSequence))
            hold on;
            idx = max([1 timeStep-100]);
            plot(eyeMovementSequence(idx:timeStep,1), -eyeMovementSequence(idx:timeStep,2), 'w-', 'LineWidth', 4.0);
            hold off;
        end
        axis 'image'; axis 'xy';
        set(gca, 'CLim', [0 1023], 'XTick', [], 'YTick', [], 'XLim', mosaicCols/2*[-1 1], 'YLim', mosaicRows/2*[-1 1], 'FontSize', 14);
        ylabel(sprintf('space\n(%2.0f retinal microns, %2.2f deg)', mosaicSize(2)*1e6, mosaicFOV(2)), 'FontSize', 16, 'FontWeight', 'bold');
        title(sprintf('t : %2.2f ms\nL-mosaic %s', timeAxis(timeStep)*1000, signalName), 'FontSize', 16);
        
        
        % M-cone mosaic subplot
        subplot('Position', [0.34 0.09 0.28 0.82]);
        imagesc(mosaicXaxis, mosaicYaxis, MconeMosaicResponse);
        % Superimpose the eye movement path to this point
        if (~isempty(eyeMovementSequence))
            hold on;
            idx = max([1 timeStep-100]);
            plot(eyeMovementSequence(idx:timeStep,1), -eyeMovementSequence(idx:timeStep,2), 'w-', 'LineWidth', 4.0);
            hold off;
        end
        axis 'image'; axis 'xy'
        set(gca, 'CLim', [0 1023], 'XTick', [], 'YTick', [], 'XLim', mosaicCols/2*[-1 1], 'YLim', mosaicRows/2*[-1 1], 'FontSize', 14);
        xlabel(sprintf('space\n(%2.0f retinal microns, %2.2f deg)', mosaicSize(1)*1e6, mosaicFOV(1)), 'FontSize', 16, 'FontWeight', 'bold');
        title(sprintf('t : %2.2f ms\nM-mosaic %s', timeAxis(timeStep)*1000, signalName), 'FontSize', 16);
        

        % S-cone mosaic subplot
        subplot('Position', [0.64 0.09 0.28 0.82]);
        imagesc(mosaicXaxis, mosaicYaxis, SconeMosaicResponse);
        % Superimpose the eye movement path to this point
        if (~isempty(eyeMovementSequence))
            hold on;
            idx = max([1 timeStep-100]);
            plot(eyeMovementSequence(idx:timeStep,1), -eyeMovementSequence(idx:timeStep,2), 'w-', 'LineWidth', 4.0);
            hold off;
        end
        axis 'image'; axis 'xy'
        set(gca, 'CLim', [0 1023], 'XTick', [], 'YTick', [], 'XLim', mosaicCols/2*[-1 1], 'YLim', mosaicRows/2*[-1 1], 'FontSize', 14);
        title(sprintf('t : %2.2f ms\nS-mosaic %s', timeAxis(timeStep)*1000, signalName), 'FontSize', 16);
        
        
        % Add colorbar
        originalPosition = get(gca, 'position');
        hCbar = colorbar('eastoutside', 'peer', gca, 'Ticks', colorbarTicks, 'TickLabels', sprintf('%2.0f\n', colorbarTickLabels));
        hCbar.Orientation = 'vertical'; 
        hCbar.Label.String = signalName;
        hCbar.FontSize = 14; 
        hCbar.FontName = 'Menlo'; 
        hCbar.Color = [0.2 0.2 0.2];
        % The addition changes the figure size, so undo this change
        newPosition = get(gca, 'position');
        set(gca,'position',[newPosition(1) newPosition(2) originalPosition(3) originalPosition(4)]);
        
        drawnow;
        writerObj.writeVideo(getframe(hFig));
    end
    writerObj.close();
end


function clut = bipolarLUT(entriesNum, flag)
    if (mod(entriesNum,2) == 0)
        m1 = entriesNum*0.5;
        r = (0:m1-1)'/max(m1-1,1);
        g = r;
        r = [r; ones(m1,1)];
        g = [g; flipud(g)];
        b = flipud(r);
    else
        m1 = floor(entriesNum*0.5);
        r = (0:m1-1)'/max(m1,1);
        g = r;
        r = [r; ones(m1+1,1)];
        g = [g; 1; flipud(g)];
        b = flipud(r);
    end
    if (flag)
        clut = [r g b]; 
    else
        clut = [r g b];
    end
end