function V1filterBank = generateV1FilterBank(spatialParams, mosaicParams, topLevelDirParams)

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
            
            filterWidth = mosaicParams.fieldOfViewDegs;
            V1filterBank = makeV1filters(spatialParams, filterWidth, coneLocsInDegs, xaxis, yaxis);
            zLevels = 0.05:0.1:1.0;
            zLevels = [-fliplr(zLevels) zLevels];
            
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', 2, ...
               'colsNum', 3, ...
               'heightMargin',   0.02, ...
               'widthMargin',    0.02, ...
               'leftMargin',     0.02, ...
               'rightMargin',    0.006, ...
               'bottomMargin',   0.01, ...
               'topMargin',      0.005);
           
            hFig = figure(1); clf;
            set(hFig, 'Position', [10 10 2200 1450]);
            
            subplot('Position', subplotPosVectors(1,1).v);
                imagesc(xaxis, yaxis, spatialModulation);
                axis 'xy';
                axis 'image';
                hold on;
                % outline mosaic extent in red
                x = mosaicParams.fieldOfViewDegs * [-0.5 0.5 0.5 -0.5 -0.5];
                y = mosaicParams.fieldOfViewDegs * [-0.5 -0.5 0.5 0.5 -0.5];
                plot(x,y, 'r-', 'LineWidth', 1.5);
                hold off
                set(gca, 'XLim', spatialParams.fieldOfViewDegs/2*[-1 1], 'YLim', spatialParams.fieldOfViewDegs/2*[-1 1]);
                xlabel('degrees');
                title('stimulus and cone mosaic');
            
             subplot('Position', subplotPosVectors(2,1).v);
                imagesc(xaxis, yaxis, spatialModulation);
                axis 'xy';
                axis 'image'
                hold on;
                plot(squeeze(coneLocsInDegs(:,1)), squeeze(coneLocsInDegs(:,2)), 'r.', 'MarkerSize', 10);
                hold off;
                set(gca, 'XLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02);
                xlabel('degrees');
                title('stimulus and cone mosaic');
                
            subplot('Position', subplotPosVectors(1,2).v);
                plot(squeeze(coneLocsInDegs(:,1)), squeeze(coneLocsInDegs(:,2)), 'r.', 'MarkerSize', 10);
                hold on
                contour(xaxis, yaxis, V1filterBank.cosPhasePoolingProfile, zLevels(zLevels>0), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
                contour(xaxis, yaxis, V1filterBank.cosPhasePoolingProfile, zLevels(zLevels<0), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
                contour(xaxis, yaxis, V1filterBank.RFprofile, [0.05 0.5], 'LineWidth', 1.5, 'Color', [0.4 0.4 0.4]);
                axis 'xy';
                axis 'image'
                hold off;
                set(gca, 'XLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YTickLabels', {});
                xlabel('degrees');
                title('cos-phase V1 filter pooling profile');
                
            subplot('Position', subplotPosVectors(1,3).v);
                plot(squeeze(coneLocsInDegs(:,1)), squeeze(coneLocsInDegs(:,2)), 'r.', 'MarkerSize', 10);
                hold on
                contour(xaxis, yaxis, V1filterBank.sinPhasePoolingProfile, zLevels(zLevels>0), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
                contour(xaxis, yaxis, V1filterBank.sinPhasePoolingProfile, zLevels(zLevels<0), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
                contour(xaxis, yaxis, V1filterBank.RFprofile, [0.05 0.5], 'LineWidth', 1.5, 'Color', [0.4 0.4 0.4]);
                axis 'xy';
                axis 'image'
                hold off;
                set(gca, 'XLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YTickLabels', {});
                xlabel('degrees');
                title('sin-phase V1 filter pooling profile');
                
            subplot('Position', subplotPosVectors(2,2).v); 
                quantizationLevels = 64;
                quantizedWeights = round(V1filterBank.cosPhasePoolingWeights*quantizationLevels);
                hold on;
                for iLevel = 1:quantizationLevels
                    idx = find(quantizedWeights == iLevel);
                    c = (0.5 + 0.5*iLevel/max(abs(quantizedWeights)))*[1 0.5 0.5];
                    plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
                    idx = find(quantizedWeights == -iLevel);
                    c = (0.5 + 0.5*iLevel/max(abs(quantizedWeights)))*[0.5 0.5 1.0];
                    plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
                end
                axis 'xy';
                axis 'image'
                hold off;
                set(gca, 'XLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YTickLabels', {});
                xlabel('degrees');
                title('cos-phase V1 filter pooling weights');
                
            subplot('Position', subplotPosVectors(2,3).v);    
                quantizationLevels = 32;
                quantizedWeights = round(V1filterBank.sinPhasePoolingWeights*quantizationLevels);
                hold on;
                for iLevel = 1:quantizationLevels
                    idx = find(quantizedWeights == iLevel);
                    c = (0.5 + 0.5*iLevel/max(abs(quantizedWeights)))*[1 0.5 0.5];
                    plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
                    idx = find(quantizedWeights == -iLevel);
                    c = (0.5 + 0.5*iLevel/max(abs(quantizedWeights)))*[0.5 0.5 1.0];
                    plot(squeeze(coneLocsInDegs(idx,1)), squeeze(coneLocsInDegs(idx,2)), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
                end
                axis 'xy';
                axis 'image'
                hold off;
                set(gca, 'XLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YLim', mosaicParams.fieldOfViewDegs/2*[-1 1]*1.02, 'YTickLabels', {});
                xlabel('degrees');
                title('sin-phase V1 filter pooling weights');
            colormap(gray(1024));
            drawnow;
        otherwise
            error('Currently generating V1 filter banks for Gabor stimuli only.');
    end
    pause;
    
end


function V1filterBank = makeV1filters(spatialParams, filterWidthDegs, coneLocsDegs, xaxisDegs, yaxisDegs)

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
end

