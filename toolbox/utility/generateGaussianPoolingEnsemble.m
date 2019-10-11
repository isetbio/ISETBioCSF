function [GaussianPoolingEnsemble, hFig] = generateGaussianPoolingEnsemble(GaussianPoolingSigmaArcMin,  mosaicParams, topLevelDirParams, visualizeSpatialScheme, thresholdParams, paramsList)

    % Load the mosaic
    coneParamsList = {topLevelDirParams, mosaicParams};
    theProgram = 't_coneCurrentEyeMovementsResponseInstances';
    rwObject = IBIOColorDetectReadWriteBasic;
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    
    coneLocsInMeters = theMosaic.coneLocsHexGridAlignedWithSerializedConeMosaicResponse();
    coneLocsDegs(:,1) = coneLocsInMeters(:,1) / theMosaic.width * theMosaic.fov(1);
    coneLocsDegs(:,2) = coneLocsInMeters(:,2) / theMosaic.height * theMosaic.fov(2);
    conesNum = size(coneLocsDegs,1);
    
    % make x-axis with a spatial resolution of 2.0 microns
    fovMicrons = max(theMosaic.fov(:)) * theMosaic.micronsPerDegree;
    xAxisMicrons = 0:2:fovMicrons/2;
    xAxisDegs = [-xAxisMicrons(end:-1:1) xAxisMicrons(2:end)]/theMosaic.micronsPerDegree;
    sigmaDegs = GaussianPoolingSigmaArcMin/60;
    
    [X,Y] = meshgrid(xAxisDegs, xAxisDegs);
    rfCoordsDegs = [X(:) Y(:)];
    
    % Find nearest cone location
    [~, idx] = pdist2(rfCoordsDegs, coneLocsDegs, 'euclidean', 'Smallest', 1);
              
    
    parfor coneIndex = 1:conesNum  
        poolingProfile = gaussianKernel(sigmaDegs,coneLocsDegs(coneIndex,:),xAxisDegs);
        poolingWeights(coneIndex,:) = poolingProfile(idx);
    end
    
    GaussianPoolingEnsemble.poolingWeights = poolingWeights;
    GaussianPoolingEnsemble.coneApertureOutlines = computeConeApertureOutlines(theMosaic);
    GaussianPoolingEnsemble.coneLocsDegs = coneLocsDegs;
    
    if (visualizeSpatialScheme)
        hFig = figure(1); clf;
        set(hFig, 'Position', [10 10 1000 1000], 'Color', [1 1 1]);
        colormap(gray);
        coneApertureOutlines = GaussianPoolingEnsemble.coneApertureOutlines;
        eccentricities = sqrt(sum(coneLocsDegs.^2,2));
        [~,coneIndex] = min(eccentricities);
        quantizationLevels = 256;
        plotQuantizedWeights(gca, squeeze(poolingWeights(coneIndex,:)), quantizationLevels, coneLocsDegs, coneApertureOutlines);
        hold on
        plotCircle(coneLocsDegs(coneIndex,1), coneLocsDegs(coneIndex,2), 2.355*sigmaDegs);
        hold off
        axis 'xy'; axis 'image';box on
        set(gca, 'FontSize', 20, 'CLim', [0 1], 'XLim', max(theMosaic.fov(:))/2*[-1 1]*1.02, 'YLim', max(theMosaic.fov(:))/2*[-1 1]*1.02, 'XTickLabels', {});
        drawnow
        
        % Save figure
        theProgram = mfilename;
        rwObject = IBIOColorDetectReadWriteBasic;
        data = 0;
        fileName = sprintf('GaussianPooling');
        paramsList{numel(paramsList)+1} = thresholdParams;
        rwObject.write(fileName, data, paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
       
       
    end
end

function plotCircle(xo,yo,diam)
    angles=1:360;
    x=diam/2*cosd(angles);
    y=diam/2*sind(angles);
    plot(xo+x,yo+y,'r-', 'LineWidth', 1.5)
end

function k=gaussianKernel(sigma,pos,x)
    [X,Y]=meshgrid(x,x);
    k=exp(-0.5*((X-pos(1))/sigma).^2).*exp(-0.5*((Y-pos(2))/sigma).^2);
end

function coneApertureOutline = computeConeApertureOutlines(theMosaic)
    
    coneLocsInMeters = theMosaic.coneLocsHexGridAlignedWithSerializedConeMosaicResponse();
    coneRadiusMicrons = 0.5*diameterForCircularApertureFromWidthForSquareAperture(theMosaic.pigment.width)*1e6;
    iTheta = (0:30:360) / 180 * pi;
    if (theMosaic.shouldCorrectAbsorptionsWithEccentricity())
        % Compute ecc-varying apertures
        dx = coneRadiusMicrons * 1e-6;
        apertureOutline.x = dx * cos(iTheta);
        apertureOutline.y = dx * sin(iTheta);

        [coneApertureOutline, ~, ~] = theMosaic.computeApertureSizes(...
            dx, [], apertureOutline, [], coneLocsInMeters(:,1), coneLocsInMeters(:,2));

        coneApertureOutline.x = coneApertureOutline.x *1e6/theMosaic.micronsPerDegree;
        coneApertureOutline.y = coneApertureOutline.y *1e6/theMosaic.micronsPerDegree;
    else
        coneRadiusDegs = coneRadiusMicrons/theMosaic.micronsPerDegree;
        coneApertureOutline.x = coneRadiusDegs * cos(iTheta);
        coneApertureOutline.y = coneRadiusDegs * sin(iTheta);
    end
end
