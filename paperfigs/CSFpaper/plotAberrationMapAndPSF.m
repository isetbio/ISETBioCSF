function plotAberrationMapAndPSF

    opticsModels = availableCustomWvfOpticsModels();
    visualizedOpticsModelIndex = 7;
    visualizedOpticsModel = opticsModels{visualizedOpticsModelIndex};
    
    targetWavelengths = 550 + [0 -70 70];
    calcPupilDiametersMM = [2 3];
    
    localDir = strrep(isetRootPath, 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    showRayTranslationMap = true;
    maxAberrationMicrons = 0.8;
    
    umPerDegree = 300;
    wavefrontSpatialSamples = 261*4+1;
    psfRange = 5;
    
    kkk = -2;
    for iW = 1:numel(targetWavelengths)
        targetWavelength = targetWavelengths(iW);
        for iP = 1:numel(calcPupilDiametersMM)
            kkk = kkk + 2;
            calcPupilDiameterMM = calcPupilDiametersMM(iP);
            % Generate the wavefront map and the optics
            [theCustomOI, Zcoeffs, theWVF] = oiWithCustomOptics(visualizedOpticsModel, ...
                    wavefrontSpatialSamples, calcPupilDiameterMM, umPerDegree, ...
                    'wavelengths', 550+[-100:5:100], ...
                    'centeringWavelength', [], ...  %  do not center PSF
                    'showTranslation',false);
    
            % Select visualize wavelength
            optics = oiGet(theCustomOI, 'optics');
            wavelengths = opticsGet(optics, 'otf wave');
            [~,iw] = min(abs(wavelengths-targetWavelength));
            targetWavelength = wavelengths(iw);
            waveOTF = opticsGet(optics,'otf data',targetWavelength);
            xSfCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
            ySfCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
            [xSfGridCyclesDeg,ySfGridCyclesDeg] = meshgrid(xSfCyclesDeg,ySfCyclesDeg);
            [xGridMinutes, yGridMinutes, wavePSF] = OtfToPsf(...
                        xSfGridCyclesDeg,ySfGridCyclesDeg,fftshift(waveOTF), ...
                        'warningInsteadOfErrorForNegativeValuedPSF', 1 ...
                 );
            xMinutes = squeeze(xGridMinutes(1,:));
            yMinutes = squeeze(yGridMinutes(:,1));  
            xx = find(abs(xMinutes) <= psfRange);
            yy = find(abs(yMinutes) <= psfRange);
            wavePSF = wavePSF(yy,xx);
            xMinutes = xMinutes(xx);
            yMinutes = yMinutes(yy);
    
            aberrationMapInMicrons = wvfGet(theWVF, 'wavefrontaberrations', targetWavelength);
            pupilSupport = wvfGet(theWVF, 'pupil spatial samples', 'mm', targetWavelength);
    
    
            export = struct(...
                'format','PDF', ...
                'name', fullfile(localDir, sprintf('AberrationMapAndPSF_%dnm%2.0fPupil.pdf', targetWavelength, calcPupilDiameterMM*10)) ...
                );
    
            panelLabels(1) = struct(...
                'text', sprintf(' %s ', char(double('A')+kkk)), ...
                'xyPos', [-3/2*0.97 3/2*0.89]);

            panelLabels(2) = struct(...
                'text', sprintf(' %s ', char(double('A')+kkk+1)), ...
                'xyPos', [xMinutes(1)*1.01 yMinutes(end)*0.915]);
    
            renderComboPlot(pupilSupport, aberrationMapInMicrons, maxAberrationMicrons, xMinutes, yMinutes, wavePSF, calcPupilDiameterMM, showRayTranslationMap, panelLabels, export);
        end
    end
end

function renderComboPlot(pupilSupport, aberrationMap, maxAberrationMap, xMinutes, yMinutes, PSF, pupilDiameter, showRayTranslationMap, panelLabels, export)
    % Render figure
    hFig = figure(1); clf;
    formatFigureForPaper(hFig, 'figureType', 'ABERRATION_MAP_PSF_COMBO');
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', 2, ...
       'heightMargin', 0.10, ...
       'widthMargin', 0.03, ...
       'leftMargin', 0.02, ...
       'rightMargin', 0.001, ...
       'bottomMargin', 0.11, ...
       'topMargin', 0.01);

    if (showRayTranslationMap)
        % Compute the gradient of the pupil phase map over some range
        % (to avoid edge artifacts), and with some stride
        stride = round(numel(pupilSupport)/100);
        [X,Y] = meshgrid(pupilSupport(1:stride:end), pupilSupport(1:stride:end));
        aberrationMapSubSampled = aberrationMap(1:stride:end, 1:stride:end);
        r = sqrt(X.^2+Y.^2);
        aberrationMapSubSampled(find(r > pupilDiameter/2)) = nan;
        [px,py] = gradient(aberrationMapSubSampled);
    end
    
    % Render the aberration map
    outline.x = 0.995*pupilDiameter/2 * cosd(0:360);
    outline.y = 0.995*pupilDiameter/2 * sind(0:360);
    aberrationMapRange = maxAberrationMap*[-1 1];
    [xx,yy] = meshgrid(pupilSupport, pupilSupport);
    in = inpolygon(xx(:), yy(:), outline.x(:), outline.y(:));
    aberrationMap(~in) = nan;
    
    ax1 = subplot('Position', subplotPosVectors(1,1).v);
    imagesc(ax1, pupilSupport, pupilSupport, aberrationMap);
    hold(ax1, 'on');
    plot(ax1,outline.x, outline.y, 'k-', 'LineWidth', 1.5);
    
    if (showRayTranslationMap)
        hold(ax1, 'on');
        % Plot the gradient vector map
        vectorScale = 2.5;
        q = quiver(ax1,X,Y, px,py, vectorScale);
        q.Color = 'k';
        q.LineWidth = 1.5;
        q.AutoScale = 'off';
        q.AutoScaleFactor = 5;
        hold(ax1, 'off');
    end
    
    xlabel(ax1,'pupil plane, x'' (mm)', 'FontWeight', 'bold');
    cmap = brewermap(1024, 'RdYlBu');
    % nan points will plot as white
    cmap(1,:) = [1 1 1];
    colormap(ax1, cmap);
    hL = colorbar('Orientation', 'Vertical', ...
        'Location', 'EastOutside');
    ylabel(hL, sprintf('optical path length difference (microns)'), ...
        'FontWeight', 'bold');
    hL.XTick = -1:0.2:1;
    
    
    set(ax1, 'CLim', aberrationMapRange, 'XLim', 1.02*3/2*[-1 1], 'YLim', 1.02*3/2*[-1 1], 'FontSize', 14);
    set(ax1, 'XTick', -(3/2):0.5:(3/2), 'YTick', -(3/2):0.5:(3/2));
    
    t = text(panelLabels(1).xyPos(1), panelLabels(1).xyPos(2), panelLabels(1).text);
    formatFigureForPaper(hFig, 'figureType', 'ABERRATION_MAP_PSF_COMBO', ...
        'theAxes', ax1,  'theText', t);
    set(ax1, 'YTickLabel', {});
    
    
    % Render the PSF
    psfPlotPos = subplotPosVectors(1,2).v;
    psfPlotPos(1) = psfPlotPos(1) + 0.03;
    ax2 = subplot('Position', psfPlotPos);
    PSF = PSF / max(abs(PSF(:)));
    contourLevels = 0:0.05:1.0;
    contourf(ax2,xMinutes, yMinutes, PSF, contourLevels);
    hold(ax2, 'on');
    midRow = floor(size(PSF,1)/2)+1;
    psfRange = xMinutes(end);
    theTicks = -5:1:5;
    plot(ax2,xMinutes, -psfRange + 1.8*psfRange * squeeze(PSF(midRow,:)), '-', 'Color', [0.2 0.6 0.9], 'LineWidth', 3.0);
    plot(ax2,yMinutes, -psfRange + 1.8*psfRange * squeeze(PSF(midRow,:)), 'b-', 'LineWidth', 1.5);
    plot(ax2,[0 0], [xMinutes(1) xMinutes(end)], 'r-', 'LineWidth', 1.0);
    plot(ax2,[yMinutes(1) yMinutes(end)], [0 0], 'r-', 'LineWidth', 1.0);
    hold(ax2, 'off');
    set(ax2, 'ZLim', [0 1], 'CLim', [0 1], 'XLim', psfRange*1.05*[-1 1], 'YLim', psfRange*1.05*[-1 1], 'FontSize', 14);
    set(ax2, 'XTick', theTicks, 'YTick', theTicks);
    xlabel(ax2,'retinal image plane, x (arc min)', 'FontWeight', 'bold');
    cmap = brewermap(1024, 'greys');
    colormap(ax2, cmap); 
    
    t = text(panelLabels(2).xyPos(1), panelLabels(2).xyPos(2), panelLabels(2).text);
    formatFigureForPaper(hFig, 'figureType', 'ABERRATION_MAP_PSF_COMBO', ...
        'theAxes', ax2,  'theText', t);
    
    set(ax2, 'YTickLabel', {});
    xtickformat(ax2, '%.0f');
    if strcmp(export.format, 'PDF')
        NicePlot.exportFigToPDF(export.name, hFig, 300);
    end
end