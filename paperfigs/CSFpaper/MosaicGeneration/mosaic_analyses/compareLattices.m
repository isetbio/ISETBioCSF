function compareLattices
    coneLocsDegsISETBio = loadISETBioMosaic();
    maxEcc = max(abs(coneLocsDegsISETBio(:)));
    coneLocsDegsBradley = loadBradleyMosaic(maxEcc);
    
    qDistISETBio = computeQuality(coneLocsDegsISETBio);
    qDistBradley = computeQuality(coneLocsDegsBradley);
    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1032 1065], 'Color', [1 1 1])
    markerSize = 3;
    subplot(2,2,1)
    plotCones(coneLocsDegsISETBio, maxEcc, 'ISETBio', [0 0 0], markerSize);
    
    subplot(2,2,2);
    plotCones(coneLocsDegsBradley, maxEcc, 'Bradley', [0 0 0], markerSize);
    
    subplot(2,2,3);
    plotQuality(qDistISETBio, 'ISETBio');

    subplot(2,2,4);
    plotQuality(qDistBradley, 'Bradley');
    NicePlot.exportFigToPDF('quality.pdf', hFig, 300);
    
    x = [-1 -0.5 0 0.5 1]*0.8*maxEcc;
    deltaDegs = 0.04;
    hFig = plotFourierAnalysis(coneLocsDegsISETBio, coneLocsDegsBradley, 'ISETBio', 'Bradley', maxEcc, x,x, deltaDegs)
    NicePlot.exportFigToPDF('Fourier.pdf', hFig, 300);
end

function hFig = plotFourierAnalysis(coneLocsA, coneLocsB, mosaicNameA, mosaicNameB, maxEcc, xx,yy, deltaDegs)

    xA = coneLocsA(:,1);
    yA = coneLocsA(:,2);
    
    xB = coneLocsB(:,1);
    yB = coneLocsB(:,2);
    
    markerSize = 6;
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1285 1112], 'Color', [1 1 1]);
    subplot(10,11, [1 2 3 4 5  12 13 14 15 16  23 24 25 26 27  34 35 36 37 38  45 46 47 48 49]);
    plotCones(coneLocsA, maxEcc, mosaicNameA, [0 0 0], markerSize); hold on;
    subsetIndex = 0;
    for l = 1:numel(yy)
            yo = -yy(l);
        for k = 1:numel(xx)
            xo = xx(k);
            idx = find(sqrt((xA-xo).^2 + (yA-yo).^2)<= deltaDegs);
            subsetIndex =  subsetIndex + 1;
            coneLocsSubsetA{subsetIndex} = coneLocsA(idx,:);
            color = [1 0 0];
            plotCones(coneLocsSubsetA{subsetIndex}, maxEcc, mosaicNameA, color, markerSize);
        end
    end
    
    subplot(10,11, [7 8 9 10 11  18 19 20 21 22  29 30 31 32 33  40 41 42 43 44  51 52 53 54 55]);
    plotCones(coneLocsB, maxEcc, mosaicNameB, [0 0 0], markerSize); hold on;
    subsetIndex = 0;
    for l = 1:numel(yy)
            yo = -yy(l);
        for k = 1:numel(xx)
            xo = xx(k);
            idx = find(sqrt((xB-xo).^2 + (yB-yo).^2)<= deltaDegs);
            subsetIndex =  subsetIndex + 1;
            coneLocsSubsetB{subsetIndex} = coneLocsB(idx,:);
            color = [1 0 0];
            plotCones(coneLocsSubsetB{subsetIndex}, maxEcc, mosaicNameB, color, markerSize);
        end
    end
    
    colormap(gray(1024));
    
    subsetIndex = 0;
    for l = 1:numel(yy)
        for k = 1:numel(xx)
            subsetIndex = subsetIndex + 1;
            base = 55 + (l-1)*11;
            coneLocsSubset = coneLocsSubsetA{subsetIndex};
            coneLocsSubset = coneLocsSubset - mean(coneLocsSubset,1);
            subplot(10,11, base + k);
            plotFourier(coneLocsSubset);
            
            base = 61 + (l-1)*11;
            coneLocsSubset = coneLocsSubsetB{subsetIndex};
            coneLocsSubset = coneLocsSubset - mean(coneLocsSubset,1);
            subplot(10,11, base + k);
            plotFourier(coneLocsSubset);
        end
    end
end


function plotFourier(coneLocsSubset)
% Generate an image of these cones by digitizing their positions in a
% grid of 2*N x 2*N pixels
    N = 64;
    maxLoc = max(abs(coneLocsSubset(:)));
    deltaLoc = maxLoc/(N-1);
    coneLocsSubset = round(coneLocsSubset/deltaLoc);
    N2 = N*8;
    mosaicImage = zeros(N2,N2);
    for k = 1:size(coneLocsSubset,1)
        mosaicImage(N2/2+coneLocsSubset(k,1), N2/2+coneLocsSubset(k,2))= 1.0;
    end

    imagesc(abs(fftshift(fft2(mosaicImage))));
    axis 'image'; axis 'ij'
    NN = 64;     
    set(gca, 'XLim', N2/2+NN*[-1 1], 'YLim', N2/2+NN*[-1 1], 'XTick', [], 'YTick', []);
    drawnow;
end



function plotQuality(qDist, mosaicName)
    qLims = [0.5 1.005]; qBins = [0.3:0.01:1.0];
    [counts,centers] = hist(qDist, qBins);
    bar(centers,counts,1)
    set(gca, 'XLim', qLims, 'YLim', [0 3500], 'XTick', [0.1:0.1:1.0], 'YTick', [0:1000:5000], 'FontSize', 14);
    axis 'square'; grid on
    xlabel('q');
    ylabel('count');
end

function plotCones(coneLocs, maxEcc, mosaicName, color, markerSize)
    plot(coneLocs(:,1), coneLocs(:,2), 'o', ...
        'MarkerFaceColor', [0.8 0.8 0.8], ...
        'MarkerEdgeColor', color, ...
        'MarkerSize', markerSize);
    set(gca, 'XLim', maxEcc*[-1 1], 'YLim', maxEcc*[-1 1], ...
        'XTick', [0.5:0.1:0.5], 'YTick', [0.5:0.1:0.5], 'FontSize', 18); 
    axis 'square';
    %title(sprintf('%s (n = %d)', mosaicName, size(coneLocs,1)));
    title(sprintf('%s', mosaicName));
end

function coneLocsDegs = loadISETBioMosaic()
    fovDegs = 0.6;
    load(sprintf('../theHexMosaic%2.2fdegs.mat',max(fovDegs)));
    whos
    coneLocsDegs = theHexMosaic.coneLocsHexGrid*1e6/theHexMosaic.micronsPerDegree;
end

function coneLocsDegs = loadBradleyMosaic(maxEcc)
    out = create_GC_mosaic(1.0);
    idx = find((abs(out.all_x(:))<maxEcc) & (abs(out.all_y(:))<maxEcc));
    coneLocsDegs = [out.all_x(idx); out.all_y(idx)]';
    close all
end

function q = computeQuality(coneLocs)
    
    triangles = delaunay(squeeze(coneLocs(:,1)),squeeze(coneLocs(:,2)));
    
    trianglesNum = size(triangles,1);
    X = coneLocs(:,1);
    Y = coneLocs(:,2);
    
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
end