function visualizeOptics

    opticsModels = availableCustomWvfOpticsModels();
    visualizedOpticsModelIndex = 7;
    visualizedOpticsModel = opticsModels{visualizedOpticsModelIndex};
    
    targetWavelengths = [450 550 650];
    
    showPupilRayMap = true;
    calcPupilDiameterMM = 3;
    umPerDegree = 300;
    wavefrontSpatialSamples = 261*2+1;
    psfRange = 3;
    otfRange = 95;
    

    [theCustomOI, Zcoeffs, theWVF] = oiWithCustomOptics(visualizedOpticsModel, ...
            wavefrontSpatialSamples, calcPupilDiameterMM, umPerDegree, ...
            'wavelengths', 550+[-100:5:100], ...
            'centeringWavelength', [], ...  %  do not center PSF
            'showTranslation',false);
    
    % Select visualize wavelength
    optics = oiGet(theCustomOI, 'optics');
    wavelengths = opticsGet(optics, 'otf wave');
    
    
    for iw = 1:numel(targetWavelengths)
        targetWavelength = targetWavelengths(iw);
        
        [~,iw] = min(abs(wavelengths-targetWavelength));
        targetWavelength = wavelengths(iw);

        xSfCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
        ySfCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
        [xSfGridCyclesDeg,ySfGridCyclesDeg] = meshgrid(xSfCyclesDeg,ySfCyclesDeg);

        waveOTF = opticsGet(optics,'otf data',targetWavelength);
        [xGridMinutes, yGridMinutes, wavePSF] = OtfToPsf(...
                xSfGridCyclesDeg,ySfGridCyclesDeg,fftshift(waveOTF), ...
                'warningInsteadOfErrorForNegativeValuedPSF', true ...
         );
        xMinutes = squeeze(xGridMinutes(1,:));
        yMinutes = squeeze(yGridMinutes(:,1));    

        hFig = figure(1); clf;
        set(hFig, 'Color', [1 1 1], 'Position', [10 10 1600 450]);
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', 4, ...
            'heightMargin', 0.10, ...
            'widthMargin', 0.03, ...
            'leftMargin', 0.02, ...
            'rightMargin', 0.001, ...
            'bottomMargin', 0.04, ...
            'topMargin', 0.02);

        ax1 = subplot('Position', subplotPosVectors(1,2).v);
        ax2 = subplot('Position', subplotPosVectors(1,1).v);

        aberrationMap = wvfGet(theWVF, 'wavefrontaberrations', targetWavelength);
        pupilFunction = wvfGet(theWVF, 'pupil function', targetWavelength);
        pupilSupport = wvfGet(theWVF, 'pupil spatial samples', 'mm', targetWavelength);

        
        plotPupilFunction(aberrationMap, pupilFunction, pupilSupport, targetWavelength, visualizedOpticsModel, ax1, ax2, ...
            'showRayMap',showPupilRayMap);

        % this is how it is actually computed
        %amp = fftshift(fft2(ifftshift(pupilFunctionData{iw}.pupilfunc)));
        %inten = (amp .* conj(amp));
        %wavePSF2 = real(inten);
       
        ax3 = subplot('Position', subplotPosVectors(1,3).v);
        plotOpticalTransferFunction(xMinutes, yMinutes, wavePSF, psfRange, ...
            'retinal space (arc min)', -10:0.5:10, sprintf('point spread function @%2.0f nm', targetWavelength), ax3);

        
        ax4 = subplot('Position', subplotPosVectors(1,4).v);
        plotOpticalTransferFunction(xSfCyclesDeg, ySfCyclesDeg, fftshift(abs(waveOTF)), otfRange, ...
            'spatial frequency (cpd)', -100:20:100, sprintf('abs(optical transfer function) @%2.0f nm', targetWavelength), ax4);


        pngFigName = sprintf('%s_%2.0f.png', visualizedOpticsModel, targetWavelength);
        NicePlot.exportFigToPNG(pngFigName, hFig, 300);
    end % iw

end

function plotPupilFunction(aberrationMap, pupilFunction, pFsupport, wavelength, opticsModel, ax1, ax2, varargin)
    p = inputParser;
    p.addParameter('showRayMap', false, @islogical);
    p.parse(varargin{:});
    
    pupilPhaseMap = angle(pupilFunction);
    cmap = hsv(1024);

    r = 1.5;
    outline.x = r * cosd(1:360);
    outline.y = r * sind(1:360);
    
    % Plot the phase map
    imagesc(ax1, pFsupport, pFsupport, pupilPhaseMap);
    [min(pupilPhaseMap(:)) max(pupilPhaseMap(:))]
    hold(ax1,'on');
    plot(ax1,outline.x, outline.y, 'k-', 'LineWidth', 2.0);
    colormap(ax1, cmap);
    axis(ax1,'image'); axis(ax1, 'xy');
    set(ax1, 'CLim', [-pi pi], 'XLim', 1.5*[-1 1], 'YLim', 1.5*[-1 1], 'FontSize', 14);
    set(ax1, 'XTick', -1.5:0.5:1.5, 'YTick', -1.5:0.5:1.5);
    set(ax1, 'YTickLabel', {});
    %phaseTicks = -pi : pi/3: pi;
    %phaseTickLabels = sprintf('%2.0f deg\n', phaseTicks/pi*180);
    %colorbar('Ticks', phaseTicks, 'TickLabels', phaseTickLabels);
    xlabel(ax1,'pupil plane (mm)', 'FontWeight', 'bold');
    box(ax1, 'on'); grid(ax1, 'on');
    title(ax1,sprintf('pupil phase map @%2.0f nm', wavelength));  
    drawnow;

    
    % Compute the gradient of the pupil phase map over some range
    % (to avoid edge artifacts), and with some stride
    pupilPosRange = 1.5; stride = 6;
    [X,Y] = meshgrid(pFsupport(1:stride:end), pFsupport(1:stride:end));
    aberrationMapSubSampled = aberrationMap(1:stride:end, 1:stride:end);
    r = sqrt(X.^2+Y.^2);
    aberrationMapSubSampled(find(r > pupilPosRange)) = nan;
    [px,py] = gradient(aberrationMapSubSampled);
    
    aberrationMapRange = 0.5*[-1 1];
    
    imagesc(ax2, pFsupport, pFsupport, aberrationMap);
    hold(ax2, 'on');
    plot(ax2,outline.x, outline.y, 'k-', 'LineWidth', 2.0);
    
    if (p.Results.showRayMap)
        hold(ax2, 'on');
        % Plot the gradient vector map
        vectorScale = 1.5;
        q = quiver(ax2,X,Y, px,py, vectorScale);
        q.Color = 'k';
        q.LineWidth = 1.0;
        q.AutoScale = 'off';
        q.AutoScaleFactor = 5;
        hold(ax2, 'off');
    end
    drawnow;
    
    cmap = brewermap(1024, 'RdYlBu');
    colormap(ax2, cmap);
    %colorbar('Orientation', 'Horizontal', 'Location', 'North');
    axis(ax2, 'image'); axis(ax2, 'xy');
    set(ax2, 'CLim', aberrationMapRange, 'XLim', 1.5*[-1 1], 'YLim', 1.5*[-1 1], 'FontSize', 14);
    set(ax2, 'XTick', -1.5:0.5:1.5, 'YTick', -1.5:0.5:1.5);
    set(ax2, 'YTickLabel', {});
    xlabel(ax2,'pupil plane (mm)', 'FontWeight', 'bold');
    box(ax2, 'on'); grid(ax2, 'on');
    title(ax2,sprintf('aberration map @%2.0f nm', wavelength));
    drawnow;    
end

    
function plotOpticalTransferFunction(xSfCyclesDeg, ySfCyclesDeg, waveOTF, otfRange, theXlabel, theTicks, theTitle, ax)
    
    xx = find(abs(xSfCyclesDeg) <= otfRange);
    yy = find(abs(ySfCyclesDeg) <= otfRange);
           
    waveOTF = waveOTF(yy,xx);
    xSfCyclesDeg = xSfCyclesDeg(xx);
    ySfCyclesDeg = ySfCyclesDeg(yy);

    waveOTF = waveOTF / max(abs(waveOTF(:)));
    contourLevels = 0:0.05:1.0;
    contourf(ax,xSfCyclesDeg, ySfCyclesDeg, waveOTF, contourLevels);
    hold(ax, 'on');
    %contourLevels = 0:0.1:1.0;
    %contour(ax,xSfCyclesDeg, ySfCyclesDeg, waveOTF, contourLevels, 'LineWidth', 1.0);
    
    
    midRow = floor(size(waveOTF,1)/2)+1;
    plot(xSfCyclesDeg, -otfRange + 1.8*otfRange * squeeze(waveOTF(midRow,:)), '-', 'Color', [0.2 0.6 0.9], 'LineWidth', 3.0);
    plot(xSfCyclesDeg, -otfRange + 1.8*otfRange * squeeze(waveOTF(midRow,:)), 'b-', 'LineWidth', 1.5);
    plot(ax,[0 0], [xSfCyclesDeg(1) xSfCyclesDeg(end)], 'r-', 'LineWidth', 1.0);
    plot(ax,[ySfCyclesDeg(1) ySfCyclesDeg(end)], [0 0], 'r-', 'LineWidth', 1.0);
    hold(ax, 'off');
    axis(ax, 'image'); axis(ax, 'xy');
    set(ax, 'ZLim', [0 1], 'CLim', [0 1], 'XLim', otfRange*1.05*[-1 1], 'YLim', otfRange*1.05*[-1 1], 'FontSize', 14);
    set(ax, 'XTick', theTicks, 'YTick', theTicks);
    set(ax, 'YTickLabel', {});
    grid(ax, 'on'); box(ax, 'on');
    cmap = brewermap(1024, 'greys');
    colormap(ax, cmap);
    xlabel(ax,theXlabel, 'FontWeight', 'bold');
    title(ax, theTitle);
end
