function compareMosaics
    
    p = getpref('IBIOColorDetect');
    IBIOColorDetectOutputBaseDir = p.outputBaseDir;
    load(fullfile(IBIOColorDetectOutputBaseDir,'[c_BanksEtAlPhotocurrentAndEyeMovements]/M_hexPacking_coneSizeUm1.5797_coneSepUmNaN_VariedConeEff_rotationDegs0_eccentricityDegs0.00_LMSdensities0.60_0.30_0.10_FOVdegs0.98x0.98_intTimeMs5_photonNoiseFrozen_osModelLinear_osTimeStepMs0.10_osNoiseFrozen_apBlur1_dark_0_0_0/matfiles/t_coneCurrentEyeMovementsResponseInstances/coneMosaic.mat'));
    
    theMosaic = theData;
    coneLocsDegs = theMosaic.coneLocsHexGrid*1e6/theMosaic.micronsPerDegree;
    maxEcc = max(abs(coneLocsDegs(:)));
    fprintf('Computing Delaunay triangularization for ISETBio mosaic\n');
    trianglesConeMosaicHex = delaunay(squeeze(coneLocsDegs(:,1)),squeeze(coneLocsDegs(:,2)));
    
    
    out = create_GC_mosaic(1.0);
    close all;
    
    idx = find((abs(out.all_x(:))<maxEcc) & (abs(out.all_y(:))<maxEcc));
    coneLocsDegsBradley = [out.all_x(idx); out.all_y(idx)]';
    fprintf('Computing Delaunay triangularization for Bradley mosaic\n');
    trianglesBradley = delaunay(squeeze(coneLocsDegsBradley(:,1)),squeeze(coneLocsDegsBradley(:,2)));
    
    fprintf('Computing quality of mosaics\n');
    [qConeMosaicHex, qDistConeMosaicHex] = computeQuality(trianglesConeMosaicHex, coneLocsDegs);
    [qBradley, qDistBradley] = computeQuality(trianglesBradley, coneLocsDegsBradley);
    
    fprintf('Plotting results\n');
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 2, ...
       'heightMargin',  0.07, ...
       'widthMargin',    0.03, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.03);

    spaceLims = 0.5*0.9*[-1 1]; spaceTicks = -0.5:0.1:0.5;
    
    qLims = [0.5 1.005]; qBins = [0.3:0.01:1.0];
    
    % Open video stream
    videoFileName = 'mosaics.mp4';
    videoOBJ = VideoWriter(videoFileName, 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 20;
    videoOBJ.Quality = 100;
    videoOBJ.open();

    hFig = figure(10); clf;
    set(hFig, 'Position', [10 10 750 780], 'Color', [1 1 1]);
    
    ax1 = subplot('Position', subplotPosVectors(1,1).v);
    plot(coneLocsDegs(:,1), coneLocsDegs(:,2), 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    hold on
    %plotTriangles(trianglesConeMosaicHex, coneLocsDegs);
    plot([0 0], [-1 1], 'r-', 'LineWidth', 1.5);
    plot([-1 1], [0 0], 'r-', 'LineWidth', 1.5);
    [xx,yy] = rectCoords(-0.3, -0.3, 0.04);
    rectPlot1 = plot(xx,yy, 'g-', 'LineWidth', 2);
    hold off
    set(gca, 'XLim', spaceLims, 'YLim', spaceLims, 'XTick', spaceTicks, 'YTick', spaceTicks, 'FontSize', 14, 'Color', [0 0 0]);
    axis 'square'
    title('ISETBio mosaic');
    
    ax2 = subplot('Position', subplotPosVectors(1,2).v);
    plot(coneLocsDegsBradley(:,1), coneLocsDegsBradley(:,2), 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    hold on
    plot([0 0], [-1 1], 'r-', 'LineWidth', 1.5);
    plot([-1 1], [0 0], 'r-', 'LineWidth', 1.5);
    %plotTriangles(trianglesBradley, coneLocsDegsBradley);
    rectPlot2 = plot(xx,yy, 'g-', 'LineWidth', 2);
    hold off
    set(gca, 'XLim', spaceLims, 'YLim', spaceLims, 'XTick', spaceTicks, 'YTick', spaceTicks, 'FontSize', 14, 'Color', [0 0 0]);
    axis 'square'
    title('Bradley mosaic');
    
    deltaDegs = 0.04;
    row = 1;
    for yo = 0.4:-0.0125:-0.4
        row = row + 1;
        if (mod(row,2) == 0)
            xx  = -0.4:0.0125:0.4;
        else
            xx = 0.4:-0.0125:-0.4;
        end
        for xo = xx
            [xx,yy] = rectCoords(xo, yo, deltaDegs);

            subplot('Position', subplotPosVectors(2,1).v);
            plotFourier(coneLocsDegs, xo, yo, deltaDegs);
            set(rectPlot1, 'XData', xx, 'YData', yy);
            set(ax1, 'XLim', xo+6*deltaDegs*[-1 1], 'YLim', yo+6*deltaDegs*[-1 1]);
            title('ISETBio mosaic');

            subplot('Position', subplotPosVectors(2,2).v);
            plotFourier(coneLocsDegsBradley, xo, yo, deltaDegs);
            set(rectPlot2, 'XData', xx, 'YData', yy);
            set(ax2, 'XLim', xo+6*deltaDegs*[-1 1], 'YLim', yo+6*deltaDegs*[-1 1]);
            title('Bradley mosaic');
            drawnow

            % Add video frame
            videoOBJ.writeVideo(getframe(hFig));
        
        end
    end
    
    % Close video stream
    videoOBJ.close();

    subplot('Position', subplotPosVectors(1,3).v);
    plotTriangles(trianglesConeMosaicHex, coneLocsDegs);
    set(gca, 'XLim', spaceLims, 'YLim', spaceLims, 'XTick', spaceTicks, 'YTick', spaceTicks, 'FontSize', 14);
    axis 'square'
    title('ISETBio mosaic');
    
    subplot('Position', subplotPosVectors(2,3).v);
    plotTriangles(trianglesBradley, coneLocsDegsBradley);
    set(gca, 'XLim', spaceLims, 'YLim', spaceLims, 'XTick', spaceTicks, 'YTick', spaceTicks, 'FontSize', 14);
    axis 'square'
    title('Bradley mosaic');
    
    subplot('Position', subplotPosVectors(1,4).v);
    [counts,centers] = hist(qDistConeMosaicHex, qBins);
    bar(centers,counts,1)
    set(gca, 'XLim', qLims, 'YLim', [0 5500], 'XTick', [0.1:0.1:1.0], 'YTick', [0:1000:5000], 'FontSize', 14);
    axis 'square'; grid on
    title('ISETBio mosaic');
    
    subplot('Position', subplotPosVectors(2,4).v);
    [counts,centers] = hist(qDistBradley, qBins);
    bar(centers,counts,1)
    set(gca, 'XLim', qLims, 'YLim', [0 5500], 'XTick', [0.1:0.1:1.0], 'YTick', [0:1000:5000], 'FontSize', 14);
    axis 'square'; grid on
    title('Bradley mosaic');
    
    drawnow
end

function [xCoords, yCoords] = rectCoords(xo, yo, delta)
    i = 0:36;
    xCoords = xo+delta*cos(i*10/180*pi);
    yCoords = yo+delta*sin(i*10/180*pi); 
end

function plotFourier(coneLocs, xo, yo, deltaDegs)

    % Select subset of cones within deltaDegs * 1.5 radius of (xo,yo)
    x = coneLocs(:,1);
    y = coneLocs(:,2);
    idx = find(sqrt((x-xo).^2 + (y-yo).^2)<= deltaDegs);
    coneLocsSubset = coneLocs(idx,:);
    fprintf('Analyzed cones: %d\n', size(coneLocsSubset,1));
    
    % zero-center the positions of this subjset of cones
    coneLocsSubset = coneLocsSubset - mean(coneLocsSubset,1);

    % Generate an image of these cones by digitizing their positions in a
    % grid of 2*N x 2*N pixels
    N = 1024;
    maxLoc = max(abs(coneLocsSubset(:)));
    deltaLoc = maxLoc/(N-1);
    coneLocsSubset = round(coneLocsSubset/deltaLoc);
    N2 = 4096*2;
    mosaicImage = zeros(N2,N2);
    for k = 1:size(coneLocsSubset,1)
        mosaicImage(N2/2+coneLocsSubset(k,1), N2/2+coneLocsSubset(k,2))= 1.0;
    end
    
    % Generate Gaussian envelope
    xn = -(N2/2):((N2/2)-1);
    yn = xn;
    [X,Y] = meshgrid(xn,yn);
    gaussianEnvelopeSigma = N2/6;
    gaussianEnvelope = exp(-0.5*(X/gaussianEnvelopeSigma).^2) .* exp(-0.5*(Y/gaussianEnvelopeSigma).^2);
    
%     figure(20);
%     subplot(1,3,1);
%     imagesc(mosaicImage)
%     axis 'square'
%     subplot(1,3,2)
%     imagesc(gaussianEnvelope)
%     axis 'image'
%     subplot(1,3,3)
%     imagesc(mosaicImage .* gaussianEnvelope)
%     axis 'image'
%     colormap(gray(1024))
%     drawnow
   
    % Multiply mosaic image with Gaussian envelope image
    mosaicImage = mosaicImage .* gaussianEnvelope;
    imagesc(abs(fftshift(fft2(mosaicImage))));
    axis 'image'
    NN = 64;
    set(gca, 'XLim', N2/2+NN*[-1 1], 'YLim', N2/2+NN*[-1 1])
    colormap(gray(1024));

    set(gca, 'FontSize', 14, 'XTick', [], 'YTick', []);

end


function plotTriangles(triangles, nodeLocs)
    trianglesNum = size(triangles,1);
    X = nodeLocs(:,1);
    Y = nodeLocs(:,2);
    for triangleIndex = 1:trianglesNum 
        plot(X(triangles(triangleIndex,:)), Y(triangles(triangleIndex,:)), 'w-', 'LineWidth', 1.0);
        if (triangleIndex == 1)
            hold on;
        end
    end
end

function [qMean,q] = computeQuality(triangles, nodeLocs)
    trianglesNum = size(triangles,1);
    X = nodeLocs(:,1);
    Y = nodeLocs(:,2);
    
    q = zeros(1,trianglesNum);
    for triangleIndex = 1:trianglesNum
        for node = 1:3
            x(node) = X(triangles(triangleIndex,node));
            y(node) = Y(triangles(triangleIndex,node));
        end 
        aLength = sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2);
        bLength = sqrt((x(1)-x(3))^2 + (y(1)-y(3))^2);
        cLength = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
        q(triangleIndex) = (bLength+cLength-aLength)*(cLength+aLength-bLength)*(aLength+bLength-cLength)/(aLength*bLength*cLength);
    end
    qMean = mean(q);
end



