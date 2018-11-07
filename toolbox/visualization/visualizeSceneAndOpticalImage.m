function visualizeSceneAndOpticalImage(backgroundScene, modulatedScene, oiBackground, oiModulated, paramsList)

    % Wether to separate plots into different figures
    separatePlotsForImageAndProfile = true;
    
    sceneLMS = sceneGet(modulatedScene, 'lms');
    sceneXYZ = sceneGet(modulatedScene, 'xyz');
    oiLMS = oiGet(oiModulated, 'lms');
    oiXYZ = oiGet(oiModulated, 'xyz');
    [sceneSRGB, sceneLRGB, maxY] = xyz2srgb(sceneXYZ);
    [oiSRGB, oiLRGB, maxY] = xyz2srgb(oiXYZ);
    
    sceneSpatialSupport = sceneGet(modulatedScene, 'spatial support');
    sceneFOV = sceneGet(modulatedScene, 'horizontalFOV');
    [xSupport, ySupport] = getXYspatialSupports(sceneSpatialSupport, sceneFOV);
    
    oiFOV = oiGet(oiModulated, 'hfov');
    oiSpatialSupport = oiGet(oiModulated, 'spatial support');
    [oiXSupport, oiYSupport] = getXYspatialSupports(oiSpatialSupport, oiFOV);
    
    sceneLMScontrast = weberContrast(sceneLMS, sceneGet(backgroundScene, 'lms'));
    oiLMScontrast = weberContrast(oiLMS, oiGet(oiBackground, 'lms'));
    
    [~,centerRow] = min(abs(ySupport));
    sceneLMScontrast = squeeze(sceneLMScontrast(centerRow,:,:));
    [~,centerRow] = min(abs(oiYSupport));
    oiLMScontrast = squeeze(oiLMScontrast(centerRow,:,:));

    halfSupport = max(xSupport(:));
    
    if (halfSupport < 0.1)
        xTicks = -0.1:0.025:0.1;
        xtickLabels = {'-0.1', '', '-0.05', '', '0', '', '+0.05', '', '+0.1'};
    elseif (halfSupport < 0.25)
        xTicks = -0.2:0.05:0.2;
        xtickLabels = {'-0.2', '', '-0.1', '', '0', '', '+0.1', '', '+0.2'};
    elseif (halfSupport < 0.4)
        xTicks = -0.4:0.1:0.4;
        xtickLabels = {'-0.4', '', '-0.2', '', '0', '', '+0.2', '', '+0.4'};
    elseif (halfSupport < 1.0)
        xTicks = -1.0:0.25:1.0;
        xtickLabels = {'-1.0', '', '-0.5', '', '0', '', '+0.5', '', '+1.0'};
    elseif (halfSupport < 2.0)
        xTicks = -2.0:0.5:2.0;
        xtickLabels = {'-2.0', '', '-1.0', '', '0', '', '+1.0', '', '+2.0'};
    elseif (halfSupport < 4.0)
        xTicks = -4.0:1.0:4.0;
        xtickLabels = {'-4.0', '', '-2.0', '', '0', '', '+2.0', '', '+4.0'};
    else
        xTicks = -10:2.5:10;
        xtickLabels = {'-10.0', '', '-5.0', '', '0', '', '+5.0', '', '+10.0'};
    end
    
    if (separatePlotsForImageAndProfile)
        hFig1 = figure(1); clf;
        set(hFig1, 'Position', [10 10 350 350], 'Color', [1 1 1]);
        subplot('Position', [0.01 0.01 0.98 0.98]);
        image(xSupport, ySupport, sceneSRGB);
        axis 'image'
        set(gca, 'XLim', [xSupport(1) xSupport(end)], ...
                 'YLim', [ySupport(1) ySupport(end)], ...
                 'XTick', xTicks, 'YTick', -10.2:0.05:0.2, ...
                 'XTickLabel', {}, 'YTickLabel', {}, ...
                 'FontSize', 16);
        grid 'on';
        drawnow;
        
        hFig2 = figure(2); clf;
        set(hFig2, 'Position', [10 10 350 350], 'Color', [1 1 1]);
        subplot('Position', [0.01 0.01 0.98 0.98]);
        image(oiXSupport, oiYSupport, oiSRGB)
        axis 'image'
        set(gca, 'XLim', [xSupport(1) xSupport(end)], ...
                 'YLim', [ySupport(1) ySupport(end)], ...
                 'XTick', xTicks, 'YTick', xTicks, ...
                 'XTickLabel', {}, 'YTickLabel', {}, ...
                 'FontSize', 16);
        grid 'on';
        
        hFig3 = figure(3); clf;
        set(hFig3, 'Position', [10 10 600 345], 'Color', [1 1 1]);
        subplot('Position', [0.15 0.235 0.84 0.71]);
        plotContrasts(xSupport, sceneLMScontrast);
        ytickLabels = {'-1', '', '-.5', '', '0', '', '+.5', '', '+1'};
        set(gca, 'XLim', [xSupport(1) xSupport(end)], ...
            'XTick', xTicks, 'YLim', [-1 1], 'YTick', -1:0.25:1,  ...
            'YTickLabel', ytickLabels, 'XTickLabel', xtickLabels, ...
            'FontSize', 28, 'LineWidth', 1.0);
        grid 'on'; box 'on';
        xlabel('\it space (deg)', 'FontWeight', 'normal', 'FontSize', 36);
        ylabel('\it stimulus contrast', 'FontWeight', 'normal', 'FontSize', 36);
        drawnow
        
        hFig4 = figure(4); clf;
        set(hFig4, 'Position', [10 10 600 345], 'Color', [1 1 1]);
        subplot('Position', [0.15 0.25 0.84 0.71]);
        plotContrasts(oiXSupport, oiLMScontrast);
        set(gca, 'XLim', [xSupport(1) xSupport(end)], ...
            'XTick', xTicks, 'YLim', [-1 1], 'YTick', -1:0.25:1,  ...
            'YTickLabel', ytickLabels, 'XTickLabel', xtickLabels, ...
            'FontSize', 28, 'LineWidth', 1.0);
        grid 'on'; box 'on';
        xlabel('\it space (deg)', 'FontWeight', 'normal', 'FontSize', 36);
        ylabel('\it retinal contrast', 'FontWeight', 'normal', 'FontSize', 36);
        drawnow
    else
        hFig = figure(1); clf;
        set(hFig, 'Position', [10 10 730 500], 'Color', [1 1 1]);
        subplot('Position', [0.08 0.335 0.44 0.68]);
        image(xSupport, ySupport, sceneSRGB);
        axis 'image'
        set(gca, 'XLim', [xSupport(1) xSupport(end)], ...
                 'YLim', [ySupport(1) ySupport(end)], ...
                 'XTick', xTicks, 'YTick', -xTicks, ...
                 'XTickLabel', {}, 'YTickLabel', {}, ...
                 'FontSize', 16);
        grid 'on';

        subplot('Position', [0.55 0.335 0.44 0.68]);
        image(oiXSupport, oiYSupport, oiSRGB)
        axis 'image'
        set(gca, 'XLim', [xSupport(1) xSupport(end)], ...
                 'YLim', [ySupport(1) ySupport(end)], ...
                 'XTick', xTicks, 'YTick', xTicks, ...
                 'XTickLabel', {}, 'YTickLabel', {}, ...
                 'FontSize', 16);
        grid 'on';

        subplot('Position', [0.08 0.095 0.44 0.24]);
        plotContrasts(xSupport, sceneLMScontrast);
        ytickLabels = {'-1.0', '', '-0.5', '', '0', '', '+0.5', '', '+1.0'};
        set(gca, 'XLim', [xSupport(1) xSupport(end)], ...
            'XTick', xTicks, 'YLim', [-1 1], 'YTick', -1:0.25:1,  ...
            'YTickLabel', ytickLabels, 'XTickLabel', xtickLabels, ...
            'FontSize', 16);
        grid 'on'; box 'on';
        xlabel('space (deg)', 'FontWeight', 'bold');
        ylabel('contrast', 'FontWeight', 'bold');

        subplot('Position', [0.55 0.095 0.44 0.24]);
        plotContrasts(oiXSupport, oiLMScontrast);
        set(gca, 'XLim', [xSupport(1) xSupport(end)], ...
            'XTick', xTicks, 'YLim', [-1 1], ...
            'YTick', -1:0.25:1, 'XTickLabel', xtickLabels, 'YTickLabel', {}, ...
            'FontSize', 16);
        grid 'on'; box 'on';
        xlabel('space (deg)', 'FontWeight', 'bold');

        drawnow;
    end
    
    
    if (~isempty(paramsList))
        % Export to PDF
        theProgram = mfilename;
        rwObject = IBIOColorDetectReadWriteBasic;
        data = 0;
        if (separatePlotsForImageAndProfile)
            fileNameForPDF = sprintf('SceneRGB');
            rwObject.write(fileNameForPDF, data, paramsList, theProgram, ...
               'type', 'NicePlotExportPNG', 'FigureHandle', hFig1, 'FigureType', 'png');
           
            fileNameForPDF = sprintf('SceneProfile');
            rwObject.write(fileNameForPDF, data, paramsList, theProgram, ...
               'type', 'NicePlotExportPDF', 'FigureHandle', hFig3, 'FigureType', 'pdf');
           
            fileNameForPDF = sprintf('OpticalImageRGB');
            rwObject.write(fileNameForPDF, data, paramsList, theProgram, ...
               'type', 'NicePlotExportPNG', 'FigureHandle', hFig2, 'FigureType', 'png');
           
            fileNameForPDF = sprintf('OpticalImageProfile');
            rwObject.write(fileNameForPDF, data, paramsList, theProgram, ...
               'type', 'NicePlotExportPDF', 'FigureHandle', hFig4, 'FigureType', 'pdf');
        else
            fileNameForPDF = sprintf('SceneRGBandOpticalImageRGB');
            rwObject.write(fileNameForPDF, data, paramsList, theProgram, ...
               'type', 'NicePlotExportPNG', 'FigureHandle', hFig, 'FigureType', 'png');
        end
    end
end

function plotContrasts(x, contrasts)
    hold on
    plotMode = 1;
    if (plotMode == 1)
        plot(x, contrasts(:,1), 'k-', 'lineWidth', 4.0);
        plot(x, contrasts(:,2), 'k-', 'lineWidth', 4.0);
        plot(x, contrasts(:,3), 'k-', 'lineWidth', 4.0);
        plot(x, contrasts(:,1), 'r-', 'lineWidth', 2.0);
        plot(x, contrasts(:,2), 'g-', 'lineWidth', 2.0);
        plot(x, contrasts(:,3), 'c-', 'lineWidth', 2.0);
    else
        p1 = plot(x, contrasts(:,1), 'r:', 'lineWidth', 3);
        p2 = plot(x, contrasts(:,2), 'g.-', 'lineWidth', 3);
        p3 = plot(x, contrasts(:,3), 'b-', 'lineWidth', 1.5);
        p1.Color(4) = 0.85;
        p2.Color(4) = 0.85;
        p3.Color(4) = 0.85;
    end
end

function [xSupport, ySupport] = getXYspatialSupports(sceneSpatialSupport, sceneFOV)
    xSupport = squeeze(sceneSpatialSupport(1,:,1));
    xSupport = xSupport / max(abs(xSupport(:))) * sceneFOV/2;
    ySupport = squeeze(sceneSpatialSupport(:,1,2));
    ySupport = ySupport / max(abs(ySupport(:))) * sceneFOV/2;
end

function contrasts = weberContrast(stimulusExcitations, backgroundExcitations)
    backgroundL = mean(mean(squeeze(backgroundExcitations(:,:,1))));
    backgroundM = mean(mean(squeeze(backgroundExcitations(:,:,2))));
    backgroundS = mean(mean(squeeze(backgroundExcitations(:,:,3))));
    contrasts(:,:,1) = (stimulusExcitations(:,:,1) - backgroundL)/backgroundL;
    contrasts(:,:,2) = (stimulusExcitations(:,:,2) - backgroundM)/backgroundM;
    contrasts(:,:,3) = (stimulusExcitations(:,:,3) - backgroundS)/backgroundS;
end


