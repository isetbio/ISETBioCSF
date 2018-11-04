function plotConeApertureMTF
    
    localExportsDir = strrep(isetRootPath, 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    localResourcesDir = strrep(isetRootPath, 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/resources');
    
    load(fullfile(localResourcesDir,'BanksMosaicApertureKernel.mat'));
    spatialSupport = oiRes*(1:size(apertureKernal,1));
    BanksMosaicAperture = struct(...
        'kernel', apertureKernal, ...
    	'spatialSupport', spatialSupport-mean(spatialSupport));
    
    load(fullfile(localResourcesDir,'HexMosaicApertureKernel.mat'));
    spatialSupport = oiRes*(1:size(apertureKernal,1));
    HexMosaicAperture = struct(...
        'kernel', apertureKernal, ...
    	'spatialSupport', spatialSupport-mean(spatialSupport));
    
    export = struct(...
        'format','PDF', ...
        'name', fullfile(localExportsDir, 'ConeApertureMTFs.pdf') ...
        );
    renderPlots(BanksMosaicAperture, HexMosaicAperture, export);
end

function renderPlots(BanksMosaicAperture, HexMosaicAperture, export)

    BanksMosaicMTF = computeApertureMTF(BanksMosaicAperture);
    HexMosaicMTF = computeApertureMTF(HexMosaicAperture);
    
    % Render figure
    hFig = figure(1); clf;
    formatFigureForPaper(hFig, 'figureType', 'CONE_APERTURE');
    
    targetSF = 60;
    [~,idx] = min(abs(BanksMosaicMTF.support1D-targetSF));
    fprintf('Banks mosaic MTF at %2.1f c/deg: %2.2f\n', BanksMosaicMTF.support1D(idx), BanksMosaicMTF.mtf1D(idx));
    
    [~,idx] = min(abs(HexMosaicMTF.support1D-targetSF));
    fprintf('Hex mosaic MTF at %2.1f c/deg: %2.2f\n', HexMosaicMTF.support1D(idx), HexMosaicMTF.mtf1D(idx));
    
    subplot('Position', [0.14 0.13 0.84 0.87]);
    plot(BanksMosaicMTF.support1D, BanksMosaicMTF.mtf1D,  'b-', 'LineWidth', 1.5);
    hold on;
    plot(HexMosaicMTF.support1D, HexMosaicMTF.mtf1D,  'r-', 'LineWidth', 1.5);
    hold off
    sfTicks = [1 3 10 30 100 300];
    yTicks = 0:0.2:1.0;
    set(gca, 'XScale', 'Log', 'XLim', [3 600], 'YLim', [0 1.01],...
         'XTick', sfTicks, 'YTick', yTicks, 'LineWidth', 0.75);
    set(gca, 'TickLength',[0.02, 0.02]);
    
    
    
    grid on; box on;
    axis 'square';
    hL = legend({'Banks ''87', 'ecc-varying'}, 'Location', 'NorthEast');
    
    t = text(3.3, 0.93, ' C ');
    formatFigureForPaper(hFig, 'figureType', 'CONE_APERTURE', 'theAxes', gca, 'theLegend', hL, 'theText', t);
    
    xlabel('\it spatial frequency (c/deg)', 'FontWeight', 'normal', 'FontSize', 28);
    ylabel('\it cone aperture MTF', 'FontWeight', 'normal', 'FontSize', 28);
    
    
    plotAperture(0.15, 0.21, BanksMosaicAperture, 3.0/2, [0 0 1]);
    plotAperture(0.37, 0.21, HexMosaicAperture, 1.5797/2, [1 0 0]);
    
    

    colormap(gray(1024));
    
    if strcmp(export.format, 'PDF')
        NicePlot.exportFigToPDF(export.name, hFig, 300);
    end
    
    if strcmp(export.format, 'PNG')
        NicePlot.exportFigToPNG(export.name, hFig, 300);
    end
    
end

function plotAperture(xo,yo, aperture, radius, lineColor)
    axes('Position', [xo yo 0.2 0.2]);
    XX = (aperture.spatialSupport(1)*1e6):0.0005:(aperture.spatialSupport(end)*1e6);
    [XX2,YY2] = meshgrid(XX,XX);
    ZZ = interp2(aperture.spatialSupport*1e6, aperture.spatialSupport*1e6, aperture.kernel, XX2,YY2, 'linear*');
    ZZ = ZZ/max(ZZ(:));
    RR = sqrt(XX2.^2+YY2.^2);
    idx = find(RR > radius);
    ZZ(idx) = 0.5;
    imagesc(XX, XX, ZZ);
    hold on;
    x = radius*0.97*cosd(0:360);
    y = radius*0.97*sind(0:360);
    plot(x,y,'-', 'Color', lineColor, 'LineWidth', 2)
    set(gca, 'CLim', [0 1], 'XLim', [-2 2], 'YLim', [-2 2], 'Color', [0.5 0.5 0.5], 'XTick', [], 'YTick', [], ...
        'XTickLabel', {}, 'YTickLabel', {}, 'FontSize', 14);
    xlabel(sprintf('%2.1f \\mum', radius*2))
    grid 'off'
    axis 'square';
end

function mosaicApertureMTF = computeApertureMTF(aperture)
    pad = 512;
    noPad = (size(aperture.kernel,1)-1)/2;
    padRange = pad+1+(-noPad:noPad);
    paddedKernel = zeros(2*pad,2*pad);
    paddedKernel(padRange, padRange) = aperture.kernel;
    
    aperture.spatialSupportDegrees = aperture.spatialSupport * 1e6/300;
    MTF = abs(fftshift(fft2(paddedKernel)));
    MTF = MTF / max(MTF(:));
    SFnyquist = 1/(2*(aperture.spatialSupportDegrees(2)-aperture.spatialSupportDegrees(1)));
    deltaSF = SFnyquist/pad;
    MTFsupport = (1:size(MTF,1))*deltaSF;
    MTFsupport = MTFsupport-mean(MTFsupport);
    idx = find(abs(MTFsupport) < 600);
    mosaicApertureMTF.support = MTFsupport(idx);
    mosaicApertureMTF.mtf = MTF(idx,idx);
    idx = find(mosaicApertureMTF.support>=0);
    mosaicApertureMTF.support1D = mosaicApertureMTF.support(idx);
    mosaicApertureMTF.mtf1D = squeeze(mosaicApertureMTF.mtf(idx(1), idx));
end
