function hFig = visualizeSpatialPoolingScheme(xaxis, yaxis, spatialModulation, ...
            spatialPoolingKernelParams, spatialPoolingFilter, coneLocsInDegs, mosaicFOVDegs, stimulusFOVDegs, coneApertureOutline)
        
        
    zLevels = [0.025:0.05:1.0];
    zLevels = [-fliplr(zLevels) zLevels];
    quantizationLevels = 1024;
    
    if (strcmp(spatialPoolingKernelParams.type, 'GaussianRF'))
        
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', 1, ...
           'heightMargin',   0.00, ...
           'widthMargin',    0.05, ...
           'leftMargin',     0.01, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.02, ...
           'topMargin',      0.005);
       
        hFig = figure(1); clf;
        set(hFig, 'Position', [10 10 850 740], 'Color', [1 1 1]);
    
        subplot('Position', subplotPosVectors(1,1).v);
        title('spatial pooling weights (cone mosaic view)');
        imagesc(xaxis, yaxis, 0.5 + 0.3*spatialModulation);
        hold on;
        plotQuantizedWeights(gca, spatialPoolingFilter.poolingWeights, quantizationLevels, coneLocsInDegs, coneApertureOutline);
        hold off;
        axis 'xy'; axis 'image'
        set(gca, 'CLim', [0 1], 'XLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'YLim', max([mosaicFOVDegs stimulusFOVDegs])/2*[-1 1]*1.02, 'XTickLabels', {});
        set(gca, 'FontSize', 20);  
        
    elseif ((strcmp(spatialPoolingKernelParams.type, 'V1QuadraturePair')) || ...
            (strcmp(spatialPoolingKernelParams.type, 'V1CosUnit')) )
        
        cMapForWeights = brewermap(256, '*RdBu');
        cMapForStim = brewermap(1024, '*Greys');
        
        if iscell(spatialPoolingFilter)
            hFig = visualizeEnsembleSpatialPoolingScheme(xaxis, yaxis, spatialModulation, ...
                spatialPoolingKernelParams, spatialPoolingFilter, ...
                coneLocsInDegs, mosaicFOVDegs, stimulusFOVDegs, coneApertureOutline, cMapForWeights);
            return;
        end
        
        if (strcmp(spatialPoolingKernelParams.type, 'V1CosUnit'))
            colsNum = 3;
            quadratureIndices = 1;
        else
            colsNum = 4;
            quadratureIndices = 2;
        end
        
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', colsNum, ...
           'heightMargin',   0.04, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.05, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.04, ...
           'topMargin',      0.01);
        
        hFig = figure(1); clf;
        formatFigureForPaper(hFig, ...
            'figureType', 'CONE_SPATIAL_POOLING');
        
        % Stimulus modulation
        tickIncrementDegs = 0.1;
        ax = subplot('Position', subplotPosVectors(1,1).v);
        imagesc(xaxis, yaxis, 0.5 + 0.3*spatialModulation);
        hold on;
        
        xRange(1) = -mosaicFOVDegs/2;
        xRange(2) = mosaicFOVDegs/2;
        yRange = xRange;
        
        %plot([xRange(1) xRange(end)], [0 0 ], 'k-', 'LineWidth', 1.0);
        %plot([0 0],[yRange(1) yRange(end)], 'k-', 'LineWidth', 1.0);
        hold off;
        formatFigureForPaper(hFig, ...
             'figureType', 'CONE_SPATIAL_POOLING', ...
             'theAxes', ax, ...
             'theFigureTitle', 'stimulus modulation' ...
        );

        midPoint = length(cMapForStim)/2
        bkgndColor = cMapForStim(midPoint,:);
        set(ax, 'CLim', [0 1], 'XLim', [xRange(1) xRange(end)], 'YLim', [yRange(1) yRange(end)]);
        set(ax, 'Color', bkgndColor, 'XTick', -0.2:0.1:0.2, 'YTick', -0.2:tickIncrementDegs:0.2);
        grid(ax, 'off')
        colormap(ax, cMapForStim);
        ylabel(ax, 'degrees');
        
        % Spatial RF envelope
        %envelopeLevels = [0.05:0.05:1];
        %contour(xaxis, yaxis, spatialPoolingFilter.RFprofile, envelopeLevels, 'Color', 'g', 'LineWidth', 1.5, 'LineStyle', '-');
   
        % The cone mosaic
        if iscell(spatialPoolingFilter)
           % ensemble of V1-receptive field like filters
           theSpatialPoolingFilter = spatialPoolingFilter{1};
        else
           theSpatialPoolingFilter = spatialPoolingFilter;
        end
            
        ax = subplot('Position', subplotPosVectors(1,2).v);
        hold on
        plotQuantizedWeights(gca, 0*theSpatialPoolingFilter.cosPhasePoolingWeights, quantizationLevels, coneLocsInDegs, coneApertureOutline);
        plot(ax,[xRange(1) xRange(end)], [0 0 ], 'k-', 'LineWidth', 1.0);
        plot(ax,[0 0],[yRange(1) yRange(end)], 'k-', 'LineWidth', 1.0);
        hold off;
        formatFigureForPaper(hFig, ...
             'figureType', 'CONE_SPATIAL_POOLING', ...
             'theAxes', ax, ...
             'theFigureTitle', 'cone mosaic' ...
        );
        set(gca, 'CLim', [0 1], 'XLim', [xRange(1) xRange(end)], 'YLim', [yRange(1) yRange(end)]);
        set(gca, 'XTick', -0.2:0.1:0.2, 'YTick', -0.2:tickIncrementDegs:0.2);
        set(gca, 'YTickLabel', {});
        colormap(gca, cMapForStim);
        %xlabel(ax, 'degrees');
        %ylabel(ax, 'degrees');
    
        if iscell(spatialPoolingFilter)
            for quadratureIndex = 1:2
                ax = subplot('Position', subplotPosVectors(1, quadratureIndex+2).v);
                imagesc(ax, xaxis, yaxis, 0.5 + 0.0*spatialModulation);
            end
        end
        
        % which bandwidth index to display
        displayedBandwidthIndex = 3;
        % which orientation index to display
        displayedOrientationIndex = 2;
                
        for unitIndex = 1:numel(spatialPoolingFilter)
            fprintf('displaying weights for unit %d of %d\n', unitIndex, numel(spatialPoolingFilter));
            
            if iscell(spatialPoolingFilter)
                % ensemble of V1-receptive field like filters
                theSpatialPoolingFilter = spatialPoolingFilter{unitIndex};
            else
                theSpatialPoolingFilter = spatialPoolingFilter;
            end
        
            if (isstruct(spatialPoolingFilter)) || ...
               ( (iscell(spatialPoolingFilter)) && ...
                 (all(theSpatialPoolingFilter.rowColPosition == [0 0])) && ...
                 ((theSpatialPoolingFilter.bandwidthIndex   == displayedBandwidthIndex)   && ...
                  (theSpatialPoolingFilter.orientationIndex == displayedOrientationIndex)) ...
               )
            
                fprintf('displaying unit %d of %d\n', unitIndex, numel(spatialPoolingFilter))
                
                % The cos/sin-weights
                maxWeight = max([max(abs(theSpatialPoolingFilter.cosPhasePoolingWeights(:))) max(abs(theSpatialPoolingFilter.sinPhasePoolingWeights(:)))]);
                for quadratureIndex = 1:quadratureIndices
                    if (quadratureIndex == 1)
                        quantizedWeights = theSpatialPoolingFilter.cosPhasePoolingWeights;
                        desiredProfile = theSpatialPoolingFilter.cosPhasePoolingProfile;
                    else
                        quantizedWeights = theSpatialPoolingFilter.sinPhasePoolingWeights;
                        desiredProfile = theSpatialPoolingFilter.sinPhasePoolingProfile;
                    end

                    ax = subplot('Position', subplotPosVectors(1, quadratureIndex+2).v);
                    hold(ax, 'on');
                    plotQuantizedWeights(gca, quantizedWeights/maxWeight, quantizationLevels, coneLocsInDegs, coneApertureOutline);
                    
                    %plot(ax,[xRange(1) xRange(end)], [0 0 ], 'k-', 'LineWidth', 1.0);
                    %plot(ax,[0 0],[yRange(1) yRange(end)], 'k-', 'LineWidth', 1.0);
                    hold(ax, 'off');
                    colormap(ax, cMapForWeights);
                end % quadrature index
                drawnow;
            end % if
        end % for unitIndex

        if iscell(spatialPoolingFilter)
            for unitIndex = 1:numel(spatialPoolingFilter)
                theSpatialPoolingFilter = spatialPoolingFilter{unitIndex};
                
                if ~((theSpatialPoolingFilter.bandwidthIndex == displayedBandwidthIndex) && ...
                      (theSpatialPoolingFilter.orientationIndex == displayedOrientationIndex))
                        continue;
                end
                
                for quadratureIndex = 1:2
                    ax = subplot('Position', subplotPosVectors(1, quadratureIndex+2).v);
                    hold(ax, 'on');
                    fprintf('Plotting RF outline for unit %d of %d\n', unitIndex, numel(spatialPoolingFilter));
                    % Spatial RF envelope
                    envelopeLevels = [0.05 0.051];
                    contour(ax, xaxis, yaxis, theSpatialPoolingFilter.RFprofile, envelopeLevels, 'Color', 'k', 'LineWidth', 1.0, 'LineStyle', '-');
                    drawnow;
                    %plot(ax, theSpatialPoolingFilter.outlineX, theSpatialPoolingFilter.outlineY, 'g-');
                end
            end
        end
        
        for quadratureIndex = 1:quadratureIndices
            if (quadratureIndex == 1)
                if (strcmp(spatialPoolingKernelParams.type, 'V1CosUnit'))
                    figTitle = 'pooling weights';
                else
                    figTitle = 'cos-phase pooling';
                end
            else
                figTitle = 'sin-phase pooling';
            end
                    
            ax = subplot('Position', subplotPosVectors(1, quadratureIndex+2).v);
            formatFigureForPaper(hFig, ...
                     'figureType', 'CONE_SPATIAL_POOLING', ...
                     'theAxes', ax, ...
                     'theFigureTitle', figTitle ...
            );
            set(ax, 'CLim', [0 1], 'XLim', [xRange(1) xRange(end)], 'YLim', [yRange(1) yRange(end)]);
            set(ax, 'XTick', -0.2:tickIncrementDegs:0.2, 'YTick', -0.2:tickIncrementDegs:0.2);
            set(ax, 'YTickLabel', {});
            %xlabel(ax, 'degrees');
        end
        
    elseif (strcmp(spatialPoolingKernelParams.type, 'V1envelope'))
       subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', 2, ...
           'heightMargin',   0.001, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.001, ...
           'topMargin',      0.001);
        
        hFig = figure(1); clf;
        set(hFig, 'Position', [10 10 1400 570], 'Color', [1 1 1]);
        
        subplot('Position', subplotPosVectors(1,1).v);
        imagesc(xaxis, yaxis, spatialModulation);
        hold on;
        % Spatial RF envelope
        envelopeLevels = [0.05:0.05:1];
        contour(xaxis, yaxis, spatialPoolingFilter.RFprofile, envelopeLevels, 'Color', 'g', 'LineWidth', 1.5, 'LineStyle', '-');
        axis 'xy'; axis 'image'
        set(gca, 'Color', [0.5 0.5 0.5], 'FontSize', 16);
        
        quantizedWeights = spatialPoolingFilter.envelopePoolingWeights;
        maxWeight = max(quantizedWeights(:));
        
        subplot('Position', subplotPosVectors(1,2).v);
        imagesc(xaxis, yaxis, spatialModulation);
        hold on;
        plotQuantizedWeights(gca, quantizedWeights/maxWeight, quantizationLevels, coneLocsInDegs, coneApertureOutline);
        hold off;
        axis 'xy'; axis 'image'
        set(gca, 'Color', [0.5 0.5 0.5], 'XTickLabel', {}, 'YTickLabel', {});
    else
        error('Unknown spatialPooling filter: %s\n', spatialPoolingKernelParams.type);
    end
      
    drawnow;
end
