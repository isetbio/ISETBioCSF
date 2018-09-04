function hFigs = visualizeEnsembleSpatialPoolingScheme(xaxis, yaxis, spatialModulation, ...
            spatialPoolingKernelParams, V1filterEnsemble, coneLocsDegs, mosaicFOVDegs, stimulusFOVDegs, coneRadiusMicrons)
        
    zLevels = [0.025:0.05:1.0];
    zLevels = [-fliplr(zLevels) zLevels];
        
    micronsPerDegree = 300;
    coneRadiusDegs = coneRadiusMicrons/micronsPerDegree;
    coneX = coneRadiusDegs * cos(2*pi*(0:30:360)/360);
    coneY = coneRadiusDegs * sin(2*pi*(0:30:360)/360);
    quantizationLevels = 256;
    
    
    unitsNum = numel(V1filterEnsemble);
    hFigs = [];
    figure();
    
    displayedHalfRows = 2;
    displayedHalfCols = 3;
    
    envelopePoolingWeights = [];
    
    bandwidthIndices = 0;
    orientationIndices = 0;
    for unitIndex = 1:unitsNum
        theSpatialPoolingFilter = V1filterEnsemble{unitIndex};
        bandwidthIndex = theSpatialPoolingFilter.bandwidthIndex;
        orientationIndex = theSpatialPoolingFilter.orientationIndex;
        if (bandwidthIndex > bandwidthIndices)
            bandwidthIndices = bandwidthIndex;
        end
         if (orientationIndex > orientationIndices)
            orientationIndices = orientationIndex;
         end
        
        ft2DIndex = theSpatialPoolingFilter.ft2Dindex;
        ft2DindexList(ft2DIndex,:) = [bandwidthIndex orientationIndex];
        
        if (isempty(envelopePoolingWeights))
            envelopePoolingWeights = zeros(100,numel(theSpatialPoolingFilter.envelopePoolingWeights));
            envelopePoolingWeights(ft2DIndex,:) = theSpatialPoolingFilter.envelopePoolingWeights;
        else
            envelopePoolingWeights(ft2DIndex,:) = envelopePoolingWeights(ft2DIndex,:) + theSpatialPoolingFilter.envelopePoolingWeights;
        end
    end
    envelopePoolingWeights = envelopePoolingWeights(1:ft2DIndex,:);
    maxWeights = max(envelopePoolingWeights(:));
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
            'rowsNum', orientationIndices, ...
            'colsNum', bandwidthIndices, ...
            'heightMargin',   0.005, ...
            'widthMargin',    0.005, ...
            'leftMargin',     0.005, ...
            'rightMargin',    0.001, ...
            'bottomMargin',   0.001, ...
            'topMargin',      0.001);
        
    hFigEnvelopes = figure(1000); clf;
    set(hFigEnvelopes, 'Position', [10 10 1000 1000], 'Color', [1 1 1]);
    
    for ft2DIndex = 1:size(envelopePoolingWeights,1)
        r = squeeze(ft2DindexList(ft2DIndex,:));
        subplot('Position', subplotPosVectors(r(2),r(1)).v);
        weights = squeeze(envelopePoolingWeights(ft2DIndex,:));
        hold on;
        plotQuantizedWeights(weights/maxWeights, quantizationLevels, coneLocsDegs, coneX, coneY);
        plot(gca,[xaxis(1) xaxis(end)], [0 0 ], 'k-', 'LineWidth', 1.0);
        plot(gca,[0 0],[yaxis(1) yaxis(end)], 'k-', 'LineWidth', 1.0);
        axis 'image'; axis 'xy';  box 'on'
        set(gca, 'Color', [0.5 0.5 0.5]);
        set(gca, 'CLim', [0 1], 'XLim', [xaxis(1) xaxis(end)], 'YLim', [yaxis(1) yaxis(end)]);
        set(gca, 'XTickLabel', {}, 'YTickLabel', {});
    end
    drawnow;
    
    for unitIndex = 1:unitsNum 
       theSpatialPoolingFilter = V1filterEnsemble{unitIndex};
       bandwidthIndex = theSpatialPoolingFilter.bandwidthIndex;
       orientationIndex = theSpatialPoolingFilter.orientationIndex;
       ft2DIndex = theSpatialPoolingFilter.ft2Dindex;
       
       if ((abs(theSpatialPoolingFilter.rowColPosition(1)) > displayedHalfRows) || ...
           (abs(theSpatialPoolingFilter.rowColPosition(2)) > displayedHalfCols))
        continue;
       end
       
       if (numel(hFigs) < ft2DIndex)
            colormap(gray(1024));
            drawnow;
            rows = theSpatialPoolingFilter.rowsNum;
            cols = theSpatialPoolingFilter.colsNum;
    
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
            'rowsNum', 2*displayedHalfRows+1, ...
            'colsNum', 2*displayedHalfCols+1, ...
            'heightMargin',   0.005, ...
            'widthMargin',    0.005, ...
            'leftMargin',     0.005, ...
            'rightMargin',    0.001, ...
            'bottomMargin',   0.001, ...
            'topMargin',      0.001);
        
            hFigs(ft2DIndex) = figure(1000+ft2DIndex); clf;
            set(hFigs(ft2DIndex), 'Position', [10+ft2DIndex*100 10+ft2DIndex*50 1650 1000], 'Color', [1 1 1]);
            %set(gcf,'renderer','opengl');
       end
       
       maxWeight = max([max(abs(theSpatialPoolingFilter.cosPhasePoolingWeights(:))) max(abs(theSpatialPoolingFilter.sinPhasePoolingWeights(:)))]);
       quantizedWeights = theSpatialPoolingFilter.cosPhasePoolingWeights;
       desiredProfile = theSpatialPoolingFilter.cosPhasePoolingProfile;
       
       row = theSpatialPoolingFilter.rowColPosition(1) + displayedHalfRows+1;
       col = theSpatialPoolingFilter.rowColPosition(2) + displayedHalfCols+1;
       subplot('Position', subplotPosVectors(displayedHalfRows*2+1+1-row,col).v);
       
       % Plot the stimulus
%        ii = 1:4:numel(xaxis);
%        jj = 1:4:numel(yaxis);
%        imagesc(xaxis(ii), yaxis(jj), 0.5 + 0.3*spatialModulation(jj,ii));
       hold on;
       plotQuantizedWeights(quantizedWeights/maxWeight, quantizationLevels, coneLocsDegs, coneX, coneY);
       plotHorizontalPoolingProfile(xaxis, min(yaxis) + 0.1*((max(yaxis)-min(yaxis))), (max(yaxis)-min(yaxis)) * 0.1, quantizedWeights, coneLocsDegs, desiredProfile, coneRadiusDegs);
                   
       plot(gca,[xaxis(1) xaxis(end)], [0 0 ], 'k-', 'LineWidth', 1.0);
       plot(gca,[0 0],[yaxis(1) yaxis(end)], 'k-', 'LineWidth', 1.0);
       axis 'image'; axis 'xy';  box 'on'
       set(gca, 'Color', [0.5 0.5 0.5]);
       set(gca, 'CLim', [0 1], 'XLim', [xaxis(1) xaxis(end)], 'YLim', [yaxis(1) yaxis(end)]);
       set(gca, 'XTickLabel', {}, 'YTickLabel', {});
    end % 
    
    colormap(gray(1024));
    drawnow;
    
    hFigs(numel(hFigs)+1) = hFigEnvelopes;
end

function plotConeLocations(coneLocsDegs, coneX, coneY, xaxis, yaxis)
    % All cone locations
    X = zeros(numel(coneX), size(coneLocsDegs,1));
    Y = X;
    for k = 1:size(coneLocsDegs,1)
        X(:,k) = coneLocsDegs(k,1)+coneX;
        Y(:,k) = coneLocsDegs(k,2)+coneY;
    end
    line(X,Y,'color','k', 'LineWidth', 0.75)
    plot([xaxis(1) xaxis(end)], [0 0 ], 'k-', 'LineWidth', 1.0);
    plot([0 0],[yaxis(1) yaxis(end)], 'k-', 'LineWidth', 1.0);
        
end

function plotQuantizedWeights(quantizedWeights, quantizationLevels, coneLocsInDegs, coneX, coneY)
            
    quantizedWeights(quantizedWeights >  1) = 1;
    quantizedWeights(quantizedWeights < -1) = -1;
    quantizedWeights = round(quantizationLevels/2 * (1+quantizedWeights));
    
    for iLevel = quantizationLevels/2:quantizationLevels
        idx = find(quantizedWeights==iLevel);
        if (~isempty(idx))
            for k = 1:numel(idx)
                c = [0.5 0.5 0.5] + (iLevel-quantizationLevels/2)/(quantizationLevels/2)*[0.5 -0.4 -0.4];
                fill(squeeze(coneLocsInDegs(idx(k),1))+coneX, squeeze(coneLocsInDegs(idx(k),2))+coneY,  c);
            end
        end
    end

    for iLevel = 0:quantizationLevels/2-1
        idx = find(quantizedWeights==iLevel);
        if (~isempty(idx))
            for k = 1:numel(idx)
                c = [0.5 0.5 0.5] + (quantizationLevels/2-iLevel)/(quantizationLevels/2)*[-0.4 -0.4 0.5];
                fill(squeeze(coneLocsInDegs(idx(k),1))+coneX, squeeze(coneLocsInDegs(idx(k),2))+coneY,  c);
            end
        end
    end
end

function plotHorizontalPoolingProfile(xaxis, y0, yA, quantizedWeights, coneLocsInDegs, desired2DProfile, coneRadiusDegs)

    xaxisIncrement = coneRadiusDegs/2;
    sampledLocations = 0:xaxisIncrement:xaxis(end);
    sampledLocationsMinus = -fliplr(sampledLocations);
    sampledLocations = [sampledLocationsMinus sampledLocations(2:end)];
    indices = zeros(1,numel(sampledLocations));
    for k = 1:numel(sampledLocations)
        [~,indices(k)] = min(abs(xaxis-sampledLocations(k)));
    end
    
    conesNum = size(coneLocsInDegs,1);
    measuredHorizontalProfile = zeros(1, numel(sampledLocations));
    for k = 1:numel(sampledLocations)
        for coneIndex = 1:conesNum
            xCoord = coneLocsInDegs(coneIndex,1);
            if (abs(xCoord-sampledLocations(k)) < coneRadiusDegs)
                measuredHorizontalProfile(k) = measuredHorizontalProfile(k) + quantizedWeights(coneIndex);
            end
        end
    end
    measuredHorizontalProfile = measuredHorizontalProfile / max(abs(measuredHorizontalProfile));

    horizontalProfile = sum(desired2DProfile,1);
    horizontalProfile = horizontalProfile/max(abs(horizontalProfile));
    
    stairs(xaxis(indices), y0 + yA * horizontalProfile(indices), 'k-', 'LineWidth', 1.5);
    stairs(xaxis(indices), y0 + yA * measuredHorizontalProfile, 'g-', 'LineWidth', 2.0);
end

