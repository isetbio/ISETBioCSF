function visualizeSceneOrOpticalImageSequence(conditionDir, sceneOrOpticalImage, sceneOrOpticalImageSequence, timeAxis, showLuminanceMap, movieFileName)
% visualizeSceneOrOpticalImageSequence(conditionDir, sceneOrOpticalImage, sceneOrOpticalImageSequence, movieFileName)
%
% Visualize a sequence of scenes or optical images and render a video of it.
%
%
%  7/12/16  npc Wrote it.

% Loop over the sequence frames and extract the needed info
for iFrame = 1:numel(sceneOrOpticalImageSequence)
    % Get an RGB rendition of the scene or the optical image
    if (strcmp(sceneOrOpticalImage, 'scene'))
        rgbFrame = sceneGet(sceneOrOpticalImageSequence{iFrame}, 'rgb image');
        if (showLuminanceMap)
            luminanceFrame = sceneGet(sceneOrOpticalImageSequence{iFrame}, 'luminance');
        end
        if (iFrame == 1)
            spatialSupport = sceneGet(sceneOrOpticalImageSequence{iFrame},'spatialsupport','cm');
            rgbImageSequence = zeros([size(rgbFrame) numel(sceneOrOpticalImageSequence)]);
            if (showLuminanceMap)
                luminanceImageSequence = zeros([size(luminanceFrame) numel(sceneOrOpticalImageSequence)]);
            end
            xAxis = squeeze(spatialSupport(1,:,1));
            yAxis = squeeze(spatialSupport(:,1,2));
        end
    elseif (strcmp(sceneOrOpticalImage, 'optical image'))
        rgbFrame = oiGet(sceneOrOpticalImageSequence{iFrame}, 'rgb image');
        if (showLuminanceMap)
            luminanceFrame = sceneGet(sceneOrOpticalImageSequence{iFrame}, 'luminance');
        end
        if (iFrame == 1)
            spatialSupport = oiGet(sceneOrOpticalImageSequence{iFrame},'spatialsupport','microns');
            rgbImageSequence = zeros([size(rgbFrame) numel(sceneOrOpticalImageSequence)]);
            if (showLuminanceMap)
                luminanceImageSequence = zeros([size(luminanceFrame) numel(sceneOrOpticalImageSequence)]);
            end
            xAxis = squeeze(spatialSupport(1,:,1));
            yAxis = squeeze(spatialSupport(:,1,2));
        end
    else
        error('First argument to visualizeSceneOrOpticalImageSequence must be one of the following strings: ''scene'' or ''optical image''\n');
    end
    
    rgbImageSequence(:,:,:,iFrame) = rgbFrame;
    if (showLuminanceMap)
        luminanceImageSequence(:,:,iFrame) = luminanceFrame;
    end
end % iFrame

% Determine ranges
rgbRange = [0 1];
if (showLuminanceMap)
    luminanceRange = mean(luminanceImageSequence(:)) * [0.9 1.1];
end

hFig = figure(10); clf;
if (showLuminanceMap)
    figSize = [10 10 950 475];
else
    figSize = [10 10 500 475];
end

set(hFig, 'Position', figSize, 'Color', [1 1 1]);
if (showLuminanceMap)
    colormap(hot(1024));
end

% Open video stream
videoDir = colorGaborDetectOutputDir(conditionDir,'videos');
videoFilename = fullfile(videoDir, sprintf('%s.m4v', movieFileName));
writerObj = VideoWriter(videoFilename, 'MPEG-4'); % H264 format
writerObj.FrameRate = 15;
writerObj.Quality = 100;
writerObj.open();


for iFrame = 1:numel(sceneOrOpticalImageSequence)
    % The RGB rendition of the scene/optical image
    if (showLuminanceMap)
        subplot('Position', [0.07 0.05 0.38 0.93]);
    else
        subplot('Position', [0.09 0.085 0.87 0.87]);
    end
    
    imagesc(xAxis, yAxis, squeeze(rgbImageSequence(:,:,:,iFrame)));
    axis 'image'
    set(gca, 'XLim', [xAxis(1) xAxis(end)], 'YLim', [yAxis(1) yAxis(end)], 'CLim', rgbRange, 'FontSize', 14);
    if (strcmp(sceneOrOpticalImage, 'scene'))
        xlabel('space (world cm)', 'FontSize', 16, 'FontWeight', 'bold');
        ylabel('space (world cm)', 'FontSize', 16, 'FontWeight', 'bold');
        title(sprintf('scene at t : %2.3f sec', timeAxis(iFrame)), 'FontSize', 16);
    else
        xlabel('space (retinal microns)', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('space (retinal microns)', 'FontSize', 14, 'FontWeight', 'bold');
        title(sprintf('optical image at t: %2.3f sec', timeAxis(iFrame)), 'FontSize', 16);
    end
    
    if (showLuminanceMap)
        % The luminance/illuminance map
        subplot('Position', [0.52 0.05 0.38 0.93]);
        imagesc(xAxis, yAxis, squeeze(luminanceImageSequence(:,:,iFrame)));
        axis 'image'
        set(gca, 'XLim', [xAxis(1) xAxis(end)], 'YLim', [yAxis(1) yAxis(end)], 'CLim', luminanceRange, 'FontSize', 14);
        if (strcmp(sceneOrOpticalImage, 'scene'))
            xlabel('space (world cm)', 'FontSize', 16, 'FontWeight', 'bold');
            title(sprintf('scene luminance at t : %2.3f sec', timeAxis(iFrame)), 'FontSize', 16);
        else
            xlabel('space (retinal microns)', 'FontSize', 14, 'FontWeight', 'bold');
            title(sprintf('retinal illuminance at t : %2.3f sec', timeAxis(iFrame)), 'FontSize', 16);
        end
        
        % Add colorbar
        originalPosition = get(gca, 'position');
        hCbar = colorbar('eastoutside', 'peer', gca); % , 'Ticks', cbarStruct.ticks, 'TickLabels', cbarStruct.tickLabels);
        hCbar.Orientation = 'vertical';
        hCbar.Label.String = 'cd/m2';
        hCbar.FontSize = 14;
        hCbar.FontName = 'Menlo';
        hCbar.Color = [0.2 0.2 0.2];
        % The addition changes the figure size, so undo this change
        newPosition = get(gca, 'position');
        set(gca,'position',[newPosition(1) newPosition(2) originalPosition(3) originalPosition(4)]);
    end
    
    % Render
    drawnow
    writerObj.writeVideo(getframe(hFig));
end % iFrame

writerObj.close();
fprintf('Movie saved in %s\n', videoFilename);

end

