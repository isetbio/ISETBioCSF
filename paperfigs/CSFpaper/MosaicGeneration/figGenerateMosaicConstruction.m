function figGenerateMosaicConstruction()

%% Initialize
ieInit; clear; close all;

% cd to wherver this script resides
[localDir,~] = fileparts(which(mfilename()));
cd(localDir)

% Set random seed to obtain replicable results
rng(1235);

params.fovDegs = [0.6 0.6]; % [1.15 1.15]; % [0.75 0.4]; % FOV in degrees ([width height], default: 0.25x0.25

makeNew = ~true;
if (makeNew)
    % Set coneMosaicHex - specific params
    params.resamplingFactor = 9;                            % 9 is good;how fine to sample the hex mosaic positions with an underlying rect mosaic                         
    params.eccBasedConeDensity = true;                      % if true, generate a mosaic where cone spacing varies with eccentricity (Curcio model)
    params.customLambda = [];                               % cone spacing in microns (only used with regular hex mosaics)
    params.rotationDegs = 0;                                % rotation of the mosaic, in degrees 0 = , 30 = (makes sense with regular hex mosaics)
    params.customInnerSegmentDiameter = [];                 % inner segment diameter, in microns (empty for isetbio default) 
    params.spatialDensity = [0 0.53 0.27 0.2];              % K/L/M/S cone densities
    params.sConeMinDistanceFactor = 2.5;                    % min distance between neighboring S-cones = f * local cone separation - to make the S-cone lattice semi-regular
    params.sConeFreeRadiusMicrons = 45;                     % radius of S-cone free retina, in microns
    params.latticeAdjustmentPositionalToleranceF =  [];      % determines cone delta movement tolerance for terminating iterative adjustment - by default this is 0.01 (here setting it lower for faster, but less acurate mosaic generation)
    params.latticeAdjustmentDelaunayToleranceF = [];        % determines position tolerance for triggering another Delaunay triangularization - by default this is 0.001 (here setting it lower for faster, but less acurate mosaic generation)
    params.maxGridAdjustmentIterations = 1500;
    params.marginF = [];

    % Generate the mosaic
    theHexMosaic = coneMosaicHex(params.resamplingFactor, ...
        'fovDegs', params.fovDegs, ...
        'eccBasedConeDensity', params.eccBasedConeDensity, ...
        'rotationDegs', params.rotationDegs, ...
        'spatialDensity', params.spatialDensity, ...
        'sConeMinDistanceFactor', params.sConeMinDistanceFactor, ...    
        'sConeFreeRadiusMicrons', params.sConeFreeRadiusMicrons, ...    
        'customLambda', params.customLambda, ...
        'customInnerSegmentDiameter', params.customInnerSegmentDiameter, ...
        'latticeAdjustmentPositionalToleranctbUeF', params.latticeAdjustmentPositionalToleranceF, ...              
        'latticeAdjustmentDelaunayToleranceF', params.latticeAdjustmentDelaunayToleranceF, ...     
        'marginF', params.marginF, ...
        'maxGridAdjustmentIterations', params.maxGridAdjustmentIterations, ...
        'saveLatticeAdjustmentProgression', true ...  
    );
    % Save the mosaic
    save(sprintf('theHexMosaic%2.2fdegs.mat',max(params.fovDegs)), 'theHexMosaic', '-v7.3');
else
    fprintf('\nLoading mosaic ... ');
    % Load the mosaic
    %load(sprintf('theHexMosaic%2.2fdegsFudge0.9.mat',max(params.fovDegs)));
    load(sprintf('theHexMosaic%2.2fdegs.mat',max(params.fovDegs)));
    fprintf('Done\n');
end

% Display mosaic info(
theHexMosaic.displayInfo();

% Show the final LMS mosaic
hFig = theHexMosaic.visualizeGrid('visualizedConeAperture', 'geometricArea', ...
    'apertureShape', 'disks', ...
    'labelConeTypes', true, ...
    'generateNewFigure', true);

cd(localDir)
NicePlot.exportFigToPDF('HexMosaic.pdf', hFig, 300);


contourLevels = 1e3 * [140 160 180 210 240 270];
    
% Show how the lattice of an ecc-based cone density hex mosaic is iteratively adjusted
hFig = theHexMosaic.plotMosaicProgression(...
    'contourLevels', contourLevels, ...
    'intermediateIterationsToDisplay', [10 100], ...
    'displayedXrangeDegs', 0.55, ...
    'displayedYrangeDegs', 0.45 ...
    );

% Export to PDF
cd(localDir)
NicePlot.exportFigToPDF('HexMosaicConstructionPartA.pdf', hFig, 300);

% Show the whole mosaic
hFig = theHexMosaic.visualizeGrid('visualizedConeAperture', 'geometricArea', ...
    'apertureShape', 'disks', ...
    'labelConeTypes', false, ...
    'coneDensityContourLevels', contourLevels, ...
    'overlayConeDensityContour', 'theoretical_and_measured', ...
    'overlayContourLabels', true, ...
    'generateNewFigure', true);
cd(localDir)
NicePlot.exportFigToPDF('HexMosaicDensity.pdf', hFig, 300);


% Show cone separation changes as the mosaic lattice converges
neigboringConesNum = 6;
hFig = visualizeConeSeparationProgression(theHexMosaic, ...
    'sampledXPositionsMicrons', [0.5 10 20 40 60 80], ...
    'neigboringConesNum', neigboringConesNum);
cd(localDir)
NicePlot.exportFigToPDF('HexMosaicConstructionPartB.pdf', hFig, 300);

end

function hFig = visualizeConeSeparationProgression(obj, varargin)

    p = inputParser;
    p.addParameter('sampledXPositionsMicrons', [0.5 5 10 20 40 80], ...
        @isnumeric);
    p.addParameter('figureLayout', [6 1], @isnumeric);
    p.addParameter('neigboringConesNum', 5, @isnumeric);
    p.parse(varargin{:});
    
    sampledXPositionsMicrons = p.Results.sampledXPositionsMicrons;
    figureLayout = p.Results.figureLayout;
    neigboringConesNum = p.Results.neigboringConesNum;
    
    hFig = figure(2);
    clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 400 500 1125]);

    rowsNum = figureLayout(1);
    colsNum = figureLayout(2);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum , ...
           'colsNum', colsNum , ...
           'heightMargin', 0.015, ...
           'widthMargin', 0.00, ...
           'leftMargin', 0.13, ...
           'rightMargin', 0.04, ...
           'bottomMargin', 0.05, ...
           'topMargin', 0.00);
       
    for idx = 1:numel(sampledXPositionsMicrons)
        trackedConePositionsMeters(idx,:) = [sampledXPositionsMicrons(idx) 0]*1e-6;  
    end
    
    [meanDist, minDist, maxDist, initialConePositions] = ...
       computeSeparationForTargetCones(obj,  trackedConePositionsMeters, neigboringConesNum);
    
    
    eccMeters = sampledXPositionsMicrons*1e-6;
    [squareSpacings, innerSegmentDiameters,~] = ...
        coneSizeReadData('eccentricity', eccMeters, 'angle', eccMeters*0);
    
    innerSegmentDiameters = diameterForCircularApertureFromWidthForSquareAperture(innerSegmentDiameters);
    theoreticalConeSpacings = squareSpacings;
    
    iterations = 1:size(meanDist,2);
    cMap = brewermap(7,'*Accent');
    cMap = cMap([2 4 3],:);
    
    markerSize = 7;
    for idx = 1:numel(sampledXPositionsMicrons)
        r = floor((idx-1)/colsNum)+1;
        c = mod(idx-1,colsNum)+1;
        subplot('Position', subplotPosVectors(r,c).v);
        color = [1 0.3 0.6]; % squeeze(cMap(1,:)); 
        edgeColor = color*0.7;
        plot(iterations, minDist(idx,:)*1e6, 'ro-',  'MarkerSize', markerSize, 'Color', edgeColor, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', color, 'LineWidth', 1.5); hold on;
        color = [0.9 0.8 0.3]; % squeeze(cMap(2,:)); 
        edgeColor = color*0.7;
        plot(iterations, meanDist(idx,:)*1e6, 'yo-', 'MarkerSize', markerSize, 'Color', edgeColor, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', color,  'LineWidth', 1.5); 
        color = [0.1 0.5 1.0]; % color = squeeze(cMap(3,:)); 
        edgeColor = color*0.7;
        plot(iterations, maxDist(idx,:)*1e6, 'bo-',  'MarkerSize', markerSize, 'Color', edgeColor, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', color,  'LineWidth', 1.5);
        %plot(iterations, innerSegmentDiameters(idx)*1e6 + zeros(size(iterations)), 'k-', 'LineWidth', 1.5);
        plot(iterations, theoreticalConeSpacings(idx)*1e6 + zeros(size(iterations)), 'k--', 'LineWidth', 1.5);
        set(gca, 'XLim', [1 size(meanDist,2)], 'XScale', 'log', 'YLim', [1.3 4.2]);
        set(gca, 'XTick', [1 3 10 30 100 300 1000], 'YTick', [1:.5:6], 'YTickLabel', sprintf(' %2.1f\n',[1:.5:6]), 'FontSize', 18, 'LineWidth', 1.0);
        grid off
        box off
        if (idx == 1)
            hL = legend(...
                {'min', ...
                sprintf('mean (n=%1.0f neighbors)', neigboringConesNum), ...
                sprintf('max (n=%1.0f neighbors)',neigboringConesNum), ...
                'Curcio ''90'}, 'Location', 'NorthWest');
            hL.NumColumns = 1;
        end
        if (r==rowsNum)
            xlabel('\it iteration no.', 'FontSize', 24);
        else
            set(gca, 'XTickLabel', {});
        end
        if (c == 1) && (r==rowsNum)
            ylabel('\it spacing (microns)','FontSize', 24)
        else
        end
        text(70, 3.7, sprintf('eccentricity: %2.0f microns', sampledXPositionsMicrons(idx)), ...
            'FontSize', 14);
    end
    
end

function [meanDist,  minDist, maxDist, initialConePositions] = computeSeparationForTargetCones(obj, trackedConePositions, neigboringConesNum)

    conePositions  = squeeze(obj.latticeAdjustmentSteps(1, :, :));
    for idx = 1:size(trackedConePositions,1)
        [~, trackedConeIndex(idx)] = min(sum((bsxfun(@minus, conePositions, trackedConePositions(idx,:))).^2, 2));
        initialConePositions(idx,:) = conePositions(trackedConeIndex(idx),:)*1e6;
    end
    
    iterationsNum = size(obj.latticeAdjustmentSteps,1);
    
    for iterIndex = 1:iterationsNum
        currentPositions = squeeze(obj.latticeAdjustmentSteps(iterIndex, :, :));
        for idx = 1:size(trackedConePositions,1)
            currentPositionsTmp = currentPositions;
            trackedConePositions(idx,:) = currentPositions(trackedConeIndex(idx),:);
            currentPositionsTmp(trackedConeIndex(idx),:) = nan;
            
            % Find the mean distances of the target cone to its neighboring cones
            nearestConeDistancesInMeters = pdist2(currentPositionsTmp, ...
            trackedConePositions(idx,:), ...
            'Euclidean', 'Smallest', neigboringConesNum);

            meanDist(idx,iterIndex) = mean(nearestConeDistancesInMeters);
            minDist(idx,iterIndex) = min(nearestConeDistancesInMeters);
            maxDist(idx,iterIndex) = max(nearestConeDistancesInMeters);
        end
        
    end
 
end




