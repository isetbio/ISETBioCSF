function V1filterBank = generateV1FilterBank(spatialParams, mosaicParams, topLevelDirParams)

    filterWidthInDegrees = 0.5*spatialParams.fieldOfViewDegs;
    plotV1FilterBank = true;
    
    % Load the mosaic
    coneParamsList = {topLevelDirParams, mosaicParams};
    theProgram = 't_coneCurrentEyeMovementsResponseInstances';
    rwObject = IBIOColorDetectReadWriteBasic;
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    
    % Get the cone locs in degrees
    if (strcmp(mosaicParams.conePacking,'rect'))
        coneLocsInMeters = theMosaic.coneLocs;
    else
        coneLocsInMeters = theMosaic.coneLocsHexGrid;
    end
    coneLocsInDegs(:,1) = coneLocsInMeters(:,1) / theMosaic.width * theMosaic.fov(1);
    coneLocsInDegs(:,2) = coneLocsInMeters(:,2) / theMosaic.height * theMosaic.fov(2);
       
    % Find the density map around each cone
    eccInMeters = sqrt(sum(coneLocsInMeters.^2, 2));
    ang = atan2(squeeze(coneLocsInMeters(:,2)), squeeze(coneLocsInMeters(:,1)))/pi*180;
    [~, ~, coneDensity] = coneSize(eccInMeters(:),ang(:));
    
    if (plotV1FilterBank)
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 3, ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.015, ...
           'topMargin',      0.005);
       
       hFig = figure(1); clf;
       set(hFig, 'Position', [10 10 2250 1450]);
    end
    
    switch(spatialParams.spatialType)
        case 'Gabor'  
            % Make the stimulus spatial modulation
            
            spatialPattern = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, 1.0));
            spatialModulation = spatialPattern-1;
            spatialModulation = spatialModulation/max(abs(spatialModulation(:)));
            
            xaxis = (0:(size(spatialModulation,2)-1))/size(spatialModulation,2) * spatialParams.fieldOfViewDegs;
            xaxis = xaxis - mean(xaxis);
            yaxis = (0:(size(spatialModulation,1)-1))/size(spatialModulation,1) * spatialParams.fieldOfViewDegs;
            yaxis = yaxis - mean(yaxis);
            
        otherwise
            error('Currently generating V1 filter banks for Gabor stimuli only.');
    end  % switch
    
    % Generate the V1 filter bank        
    V1filterBank = makeV1filters(spatialParams, filterWidthInDegrees, coneLocsInDegs, xaxis, yaxis, coneDensity);

    % Plot it
    if (plotV1FilterBank)
        zLevels = 0.05:0.2:1.0;
        zLevels = [-fliplr(zLevels) zLevels];

        subplot('Position', subplotPosVectors(1,1).v);
            imagesc(xaxis, yaxis, spatialModulation);
            axis 'xy';
            axis 'image';
            hold on;
            % V1 RF envelope in cyan
            contour(xaxis, yaxis, V1filterBank.RFprofile, zLevels(zLevels>0), 'Color', 'c', 'LineWidth', 1.0, 'LineStyle', '-');
            % outline mosaic extent in green
            x = mosaicParams.fieldOfViewDegs * [-0.5 0.5 0.5 -0.5 -0.5];
            y = mosaicParams.fieldOfViewDegs * [-0.5 -0.5 0.5 0.5 -0.5];
            plot(x,y, 'g-', 'LineWidth', 1.5);
            hold off
            set(gca, 'XLim', spatialParams.fieldOfViewDegs/2*[-1 1], 'YLim', spatialParams.fieldOfViewDegs/2*[-1 1]);
            set(gca, 'FontSize', 14, 'XTickLabels', {});
            ylabel('degrees');
            title('stimulus, cone mosaic, and V1 RF profile (stimulus view)');

         subplot('Position', subplotPosVectors(2,1).v);
            imagesc(xaxis, yaxis, spatialModulation);
            axis 'xy';
            axis 'image'
            hold on;
            plot(squeeze(coneLocsInDegs(:,1)), squeeze(coneLocsInDegs(:,2)), 'r.', 'MarkerSize', 10);
            hold off;
            set(gca, 'XLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02);
            set(gca, 'FontSize', 14, 'XTickLabels', {});
            ylabel('degrees');
            title('stimulus and cone mosaic (cone mosaic view)', 'FontSize', 16);

        subplot('Position', subplotPosVectors(1,2).v);
            imagesc(xaxis, yaxis, spatialModulation);
            axis 'xy';
            axis 'image';
            hold on;
            % outline mosaic extent in green
            x = mosaicParams.fieldOfViewDegs * [-0.5 0.5 0.5 -0.5 -0.5];
            y = mosaicParams.fieldOfViewDegs * [-0.5 -0.5 0.5 0.5 -0.5];
            plot(x,y, 'g-', 'LineWidth', 1.5);
            % V1 cos-phase filter
            contour(xaxis, yaxis, V1filterBank.cosPhasePoolingProfile, zLevels(zLevels>0), 'Color', 'r', 'LineWidth', 1.0, 'LineStyle', '-');
            contour(xaxis, yaxis, V1filterBank.cosPhasePoolingProfile, zLevels(zLevels<0), 'Color', 'b', 'LineWidth', 1.0, 'LineStyle', '-');
            plot(x,y, 'g-', 'LineWidth', 1.5);
            hold off
            set(gca, 'XLim', spatialParams.fieldOfViewDegs/2*[-1 1], 'YLim', spatialParams.fieldOfViewDegs/2*[-1 1], 'XTickLabels', {});
            set(gca, 'FontSize', 16);
            title('stimulus, cone mosaic and cos-phase V1 filter profile (stimulus view)');

        subplot('Position', subplotPosVectors(1,3).v);
            imagesc(xaxis, yaxis, spatialModulation);
            axis 'xy';
            axis 'image';
            hold on;
            % outline mosaic extent in green
            x = mosaicParams.fieldOfViewDegs * [-0.5 0.5 0.5 -0.5 -0.5];
            y = mosaicParams.fieldOfViewDegs * [-0.5 -0.5 0.5 0.5 -0.5];
            plot(x,y, 'g-', 'LineWidth', 1.5);
            % V1 sin-phase filter
            contour(xaxis, yaxis, V1filterBank.sinPhasePoolingProfile, zLevels(zLevels>0), 'Color', 'r', 'LineWidth', 1.0, 'LineStyle', '-');
            contour(xaxis, yaxis, V1filterBank.sinPhasePoolingProfile, zLevels(zLevels<0), 'Color', 'b', 'LineWidth', 1.0, 'LineStyle', '-');
            plot(x,y, 'g-', 'LineWidth', 1.5);
            hold off
            set(gca, 'XLim', spatialParams.fieldOfViewDegs/2*[-1 1], 'YLim', spatialParams.fieldOfViewDegs/2*[-1 1], 'XTickLabels', {}, 'YTickLabels', {});
            set(gca, 'FontSize', 14);
            title('stimulus, cone mosaic, and sin-phase V1 filter profile (stimulus view)');


        quantizationLevels = 1024;
        subplot('Position', subplotPosVectors(2,2).v); 
            quantizedWeights = round(V1filterBank.cosPhasePoolingWeights*quantizationLevels);
            hold on;
            for iLevel = 1:quantizationLevels
                idx = find(quantizedWeights == iLevel);
                c = (iLevel/max(abs(quantizedWeights)))*[1 0.5 0.5];
                plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
                idx = find(quantizedWeights == -iLevel);
                c = (iLevel/max(abs(quantizedWeights)))*[0.5 0.5 1.0];
                plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
            end
            axis 'xy';
            axis 'image'
            hold off;
            set(gca, 'Color', [0 0 0], 'XLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'XTickLabels', {});
            set(gca, 'FontSize', 14);
            title('cos-phase V1 filter pooling weights (cone mosaic view)');

        subplot('Position', subplotPosVectors(2,3).v);    
            quantizedWeights = round(V1filterBank.sinPhasePoolingWeights*quantizationLevels);
            hold on;
            for iLevel = 1:quantizationLevels
                idx = find(quantizedWeights == iLevel);
                c = (iLevel/max(abs(quantizedWeights)))*[1 0.5 0.5];
                plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
                idx = find(quantizedWeights == -iLevel);
                c = (iLevel/max(abs(quantizedWeights)))*[0.5 0.5 1.0];
                plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
            end
            axis 'xy';
            axis 'image'
            hold off;
            set(gca, 'Color', [0 0 0], 'XLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'XTickLabels', {}, 'YTickLabels', {});
            set(gca, 'FontSize', 14);
            title('sin-phase V1 filter pooling weights (cone mosaic view)');
        colormap(gray(1024));
        drawnow;
    end % plotFilterBank
           
end


function V1filterBank = makeV1filters(spatialParams, filterWidthDegs, coneLocsDegs, xaxisDegs, yaxisDegs, coneDensity)

    spatialParams.gaussianFWHMDegs = filterWidthDegs/2.0;
    
    % make the cos-phase filter
    spatialParams.ph = 0;
    cosPhaseFilter = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, 1.0));
    V1filterBank.cosPhasePoolingProfile = cosPhaseFilter-1;
    
    % make the sin-phase filter
    spatialParams.ph = pi/2;
    sinPhaseFilter = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, 1.0));
    V1filterBank.sinPhasePoolingProfile = sinPhaseFilter-1;

    % Compute energy envelope
    RFprofile = sqrt(V1filterBank.cosPhasePoolingProfile.^2 + V1filterBank.sinPhasePoolingProfile.^2);
    V1filterBank.RFprofile = RFprofile / max(abs(RFprofile(:)));
    
    % Find pooling weights
    [X,Y] = meshgrid(xaxisDegs, yaxisDegs);
    rfCoordsDegs = [X(:) Y(:)];
    
    [~, idx] = pdist2(rfCoordsDegs, coneLocsDegs, 'euclidean', 'Smallest', 1);
    V1filterBank.cosPhasePoolingWeights = V1filterBank.cosPhasePoolingProfile(idx);
    V1filterBank.sinPhasePoolingWeights = V1filterBank.sinPhasePoolingProfile(idx);
    
    % Adjust weights by the inverse of the coneDensity
    maxCos = max(abs(V1filterBank.cosPhasePoolingWeights(:)));
    maxSin = max(abs(V1filterBank.sinPhasePoolingWeights(:)));
    V1filterBank.cosPhasePoolingWeights = V1filterBank.cosPhasePoolingWeights ./ coneDensity;
    V1filterBank.sinPhasePoolingWeights = V1filterBank.sinPhasePoolingWeights ./ coneDensity;
    
    V1filterBank.cosPhasePoolingWeights = V1filterBank.cosPhasePoolingWeights / max(abs(V1filterBank.cosPhasePoolingWeights(:))) * maxCos;
    V1filterBank.sinPhasePoolingWeights = V1filterBank.sinPhasePoolingWeights / max(abs(V1filterBank.sinPhasePoolingWeights(:))) * maxSin;
end

