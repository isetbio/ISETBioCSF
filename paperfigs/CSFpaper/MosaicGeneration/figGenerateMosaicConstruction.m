function figGenerateMosaicConstruction()

%% Initialize
ieInit; clear; close all;

% cd to wherver this script resides
[localDir,~] = fileparts(which(mfilename()));
cd(localDir)

% Set random seed to obtain replicable results
rng(1235);

params.fovDegs = [0.6 0.6]; % [1.15 1.15]; % [0.75 0.4]; % FOV in degrees ([width height], default: 0.25x0.25

saveLatticeAdjustmentProgression = true;                % set to true, only if interested to see how the mosaic lattice is iteratively adjusted when eccBasedConeDensity is true               

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
        'saveLatticeAdjustmentProgression', saveLatticeAdjustmentProgression ...  
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

visualizeIterativeGridAdjustment(theHexMosaic);


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

hFig = theHexMosaic.visualizeGrid('visualizedConeAperture', 'geometricArea', ...
    'apertureShape', 'disks', ...
    'labelConeTypes', false, ...
    'coneDensityContourLevels', contourLevels, ...
    'overlayConeDensityContour', 'theoretical_and_measured', ...
    'overlayContourLabels', true, ...
    'generateNewFigure', true);
cd(localDir)
NicePlot.exportFigToPDF('HexMosaicDensity.pdf', hFig, 300);
    
if (saveLatticeAdjustmentProgression)
    % Show how the lattice of an ecc-based cone density hex mosaic is iteratively adjusted
    
    hFig = theHexMosaic.plotMosaicProgression(...
        'contourLevels', contourLevels, ...
        'intermediateIterationsToDisplay', [10 100], ...
        'displayedXrangeDegs', 0.55, ...
        'displayedYrangeDegs', 0.45 ...
        );

    % Export to PDF
    cd(localDir)
    NicePlot.exportFigToPDF('HexMosaicConstruction.pdf', hFig, 300);
end
end

function visualizeIterativeGridAdjustment(theHexMosaic)

    xDistances = [0.5 5 10 15 20 30 40 60 80];
    for idx = 1:numel(xDistances)
        trackedConePositions(idx,:) = [xDistances(idx) 0]*1e-6;  
    end
    
    [meanDist, minDist, maxDist, initialConePositions] = plotSeparationForTargetCones(theHexMosaic,  trackedConePositions);
    
    rowsNum = 3;
    colsNum = 3;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum , ...
           'colsNum', colsNum , ...
           'heightMargin', 0.05, ...
           'widthMargin', 0.03, ...
           'leftMargin', 0.05, ...
           'rightMargin', 0.02, ...
           'bottomMargin', 0.06, ...
           'topMargin', 0.02);
       
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1200 1100], 'Color', [1 1 1]);
    ecc = xDistances*1e-6;
    [squareSpacings,innerSegmentDiameters,~] = coneSizeReadData('eccentricity', ecc, 'angle', ecc*0);
    
    innerSegmentDiameters = diameterForCircularApertureFromWidthForSquareAperture(innerSegmentDiameters);
    theoreticalConeSpacings = squareSpacings;
    
    iterations = 1:size(meanDist,2);
    for idx = 1:numel(xDistances)
        r = floor((idx-1)/3)+1;
        c = mod(idx-1,3)+1;
        subplot('Position', subplotPosVectors(r,c).v)
        plot(iterations, meanDist(idx,:)*1e6, 'ko-', 'LineWidth', 1.5); hold on;
        plot(iterations, minDist(idx,:)*1e6, 'ro-',  'LineWidth', 1.5);
        plot(iterations, maxDist(idx,:)*1e6, 'bo-',  'LineWidth', 1.5);
        plot(iterations, innerSegmentDiameters(idx)*1e6 + zeros(size(iterations)), 'k-', 'LineWidth', 1.5);
        plot(iterations, theoreticalConeSpacings(idx)*1e6 + zeros(size(iterations)), 'k--', 'LineWidth', 1.5);
        set(gca, 'XLim', [1 size(meanDist,2)], 'XScale', 'log', 'YLim', [1.5 4.0]);
        set(gca, 'XTick', [1 3 10 30 100 300 1000], 'YTick', [1:0.5:6], 'FontSize', 18);
        grid on
        if (idx == 1)
        legend({'mean (6 neighbors)', 'min', 'max (6 neighbors)', 'inner segment diameter', 'Curcio ''90 spacing'});
        end
        if (r==3)
            xlabel('\it iteration', 'FontSize', 20);
        else
            set(gca, 'XTickLabel', {});
        end
        if (c == 1)
            ylabel('\it cone spacing (microns)','FontSize', 20)
        else
            set(gca, 'YTickLabel', {});
        end
        title(sprintf('cone eccentricity: %2.0f microns', xDistances(idx)));
    end
    
end

function [meanDist,  minDist, maxDist, initialConePositions] = plotSeparationForTargetCones(obj, trackedConePositions)

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
            
            
            % Find the mean distances of the target cone to its 6 neighboring cones
            neigboringConesNum = 6;
            nearestConeDistancesInMeters = pdist2(currentPositionsTmp, ...
            trackedConePositions(idx,:), ...
            'Euclidean', 'Smallest', neigboringConesNum);

            meanDist(idx,iterIndex) = mean(nearestConeDistancesInMeters);
            minDist(idx,iterIndex) = min(nearestConeDistancesInMeters);
            maxDist(idx,iterIndex) = max(nearestConeDistancesInMeters);
        end
        
    end
 
end


