function hFigs = visualizeEnsembleSpatialPoolingScheme(xaxis, yaxis, spatialModulation, ...
            spatialPoolingKernelParams, V1filterEnsemble, coneLocsDegs, mosaicFOVDegs, stimulusFOVDegs, coneApertureOutline, cMapForWeights)
        
    zLevels = [0.025:0.05:1.0];
    zLevels = [-fliplr(zLevels) zLevels];
        
    quantizationLevels = 256;
    
    
    unitsNum = numel(V1filterEnsemble);
    hFigs = [];
   
    
    envelopePoolingWeights = [];
    
    bandwidthIndices = 0;
    orientationIndices = 0;
    for unitIndex = 1:unitsNum
        theSpatialPoolingFilter = V1filterEnsemble{unitIndex};
        
        displayedHalfRows = max(theSpatialPoolingFilter.rowColPosition(:));
        displayedHalfCols = displayedHalfRows;
    
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
        
        envelopeProfiles(unitIndex,:,:) = theSpatialPoolingFilter.RFprofile;
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
    
    [X,Y] = meshgrid(xaxis-max(abs(coneApertureOutline.x(:)))/2, yaxis-max(abs(coneApertureOutline.y(:)))/2);
    for unitIndex = 1:unitsNum 
        if (unitIndex == 1)
            theSpatialPoolingFilter = V1filterEnsemble{unitIndex};
            quantizedWeights = theSpatialPoolingFilter.cosPhasePoolingWeights;
            maxWeight = max(quantizedWeights(:));
            plotQuantizedWeights(gca, quantizedWeights/maxWeight, quantizationLevels, coneLocsDegs, coneApertureOutline); hold on;
            colormap(gca, cMapForWeights);
        end
        
        theProfile = squeeze(envelopeProfiles(unitIndex,:,:));
        if (unitIndex == 1)
            [~,c] = contour(gca,X,Y,theProfile, 0.05*[1 1], 'k-');
        else
            [~,c] = contour(gca,X,Y,theProfile, 0.05*[1 1], 'k-', 'Color', [0.3 0.3 0.3]);
        end
        c.LineWidth = 2.0;

        
        if (ft2DIndex == unitsNum)
            plot(gca,[xaxis(1) xaxis(end)], [0 0 ], 'k-', 'LineWidth', 1.0);
            plot(gca,[0 0],[yaxis(1) yaxis(end)], 'k-', 'LineWidth', 1.0);
        end
        axis 'image'; axis 'xy';  box 'on'
        set(gca, 'Color', [1 1 1]);
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
            set(hFigs(ft2DIndex), 'Position', [10+ft2DIndex*100 10+ft2DIndex*50 1000 1000], 'Color', [1 1 1]);
            %set(gcf,'renderer','opengl');
       end
       
       
       row = theSpatialPoolingFilter.rowColPosition(1) + displayedHalfRows+1;
       col = theSpatialPoolingFilter.rowColPosition(2) + displayedHalfCols+1;
       subplot('Position', subplotPosVectors(displayedHalfRows*2+1+1-row,col).v);
       
       maxWeight = max([max(abs(theSpatialPoolingFilter.cosPhasePoolingWeights(:))) max(abs(theSpatialPoolingFilter.sinPhasePoolingWeights(:)))]);
       quantizedWeights = theSpatialPoolingFilter.cosPhasePoolingWeights;
       desiredProfile = theSpatialPoolingFilter.cosPhasePoolingProfile;
       
       plotQuantizedWeights(gca, quantizedWeights/maxWeight, quantizationLevels, coneLocsDegs, coneApertureOutline);
       hold (gca, 'on')
       plot(gca,2.0*[xaxis(1) xaxis(end)], [0 0 ], 'k-', 'LineWidth', 1.0);
       plot(gca,[0 0],2.0*[yaxis(1) yaxis(end)], 'k-', 'LineWidth', 1.0);
       hold (gca, 'off')
       colormap(gca, cMapForWeights);
       axis 'image'; axis 'xy';  box 'off'
       set(gca, 'Color', [1 1 1]);
       set(gca, 'CLim', [0 1], 'XLim', [xaxis(1) xaxis(end)], 'YLim', [yaxis(1) yaxis(end)]);
       set(gca, 'XTickLabel', {}, 'YTickLabel', {});
    end %
    drawnow;
    
    hFigs(numel(hFigs)+1) = hFigEnvelopes;
end

