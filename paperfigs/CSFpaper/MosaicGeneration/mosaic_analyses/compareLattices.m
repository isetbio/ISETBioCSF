function compareLattices
    coneLocsDegsISETBio = loadISETBioMosaic();
    maxEcc = max(abs(coneLocsDegsISETBio(:)));
    coneLocsDegsBradley = loadBradleyMosaic(maxEcc);
    
    qDistISETBio = computeQuality(coneLocsDegsISETBio);
    qDistBradley = computeQuality(coneLocsDegsBradley);
    
    figure(1); clf;
    subplot(2,2,1)
    plotCones(coneLocsDegsISETBio, maxEcc, 'ISETBio');
    
    subplot(2,2,2);
    plotCones(coneLocsDegsBradley, maxEcc, 'Bradley');
    
    subplot(2,2,3);
    plotQuality(qDistISETBio, 'ISETBio');

    subplot(2,2,4);
    plotQuality(qDistBradley, 'Bradley');
end

function plotFourier(coneLocs, xo, yo, deltaDegs)

    figure(2); clf;
    
    % Select subset of cones within deltaDegs * 1.5 radius of (xo,yo)
    x = coneLocs(:,1);
    y = coneLocs(:,2);
    
    idx = find(sqrt((x-xo).^2 + (y-yo).^2)<= deltaDegs);
    coneLocsSubset = coneLocs(idx,:);
    fprintf('Analyzed cones: %d\n', size(coneLocsSubset,1));
    
    % zero-center the positions of this subset of cones
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



function plotQuality(qDist, mosaicName)
    qLims = [0.5 1.005]; qBins = [0.3:0.01:1.0];
    [counts,centers] = hist(qDist, qBins);
    bar(centers,counts,1)
    set(gca, 'XLim', qLims, 'YLim', [0 3500], 'XTick', [0.1:0.1:1.0], 'YTick', [0:1000:5000], 'FontSize', 14);
    axis 'square'; grid on
    xlabel('q');
    ylabel('count');
end

function plotCones(coneLocs, maxEcc, mosaicName)
    plot(coneLocs(:,1), coneLocs(:,2), 'k.');
    set(gca, 'XLim', maxEcc*[-1 1], 'YLim', maxEcc*[-1 1], ...
        'XTick', [0.5:0.1:0.5], 'YTick', [0.5:0.1:0.5], 'FontSize', 18); 
    axis 'square';
    title(sprintf('%s (n = %d)', mosaicName, size(coneLocs,1)));
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