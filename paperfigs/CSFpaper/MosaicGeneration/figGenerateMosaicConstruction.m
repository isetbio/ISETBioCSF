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
% Show the mosaic with contour levels
hFig = theHexMosaic.visualizeGrid('visualizedConeAperture', 'geometricArea', ...
    'apertureShape', 'disks', ...
    'labelConeTypes', false, ...
    'coneDensityContourLevels', contourLevels, ...
    'overlayConeDensityContour', 'theoretical_and_measured', ...
    'overlayContourLabels', true, ...
    'generateNewFigure', true);
cd(localDir)
NicePlot.exportFigToPDF('HexMosaicDensity.pdf', hFig, 300);

  
hFig = visualizeLocalForceProgression(theHexMosaic);


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



% Show cone separation changes as the mosaic lattice converges
neigboringConesNum = 6;
hFig = visualizeConeSeparationProgression(theHexMosaic, ...
    'sampledXPositionsMicrons', [0.5 10 20 40 60 80], ...
    'neigboringConesNum', neigboringConesNum);
% Export to PDF
cd(localDir)
NicePlot.exportFigToPDF('HexMosaicConstructionPartC.pdf', hFig, 300);
end


function hFig = visualizeLocalForceProgression(obj, varargin)
    hFig = figure(2);
    clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 400 500 1125]);

    rowsNum = 6;
    colsNum = 1;
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum , ...
           'colsNum', colsNum , ...
           'heightMargin', 0.015, ...
           'widthMargin', 0.00, ...
           'leftMargin', 0.13, ...
           'rightMargin', 0.04, ...
           'bottomMargin', 0.05, ...
           'topMargin', 0.00);
   
    params.latticeAdjustmentPositionalToleranceF = 0.0013/8;
    params.latticeAdjustmentDelaunayToleranceF = 0.0013/8;
    params.maxGridAdjustmentIterations = Inf;
    params.saveLatticeAdjustmentProgression = true;
    params.latticeAdjustmentSteps = [];
    
    grid.lambdaMin = 2;
    grid.lambdaMid = 2;
    grid.coneSpacingFunction = @coneSpacingFunction;
    grid.domainFunction = @ellipticalDomainFunction;
    grid.center = obj.center * 1e6;
    grid.rotationAngle = obj.rotationDegs / 180 * pi;
    grid.width = obj.width * 1e6;
    grid.height = obj.height * 1e6;
    grid.radius = obj.marginF * sqrt(2) * ...
        max([grid.width / 2, grid.height / 2]);
    grid.ellipseAxes = determineEllipseAxesLength(grid.radius);
    grid.borderTolerance = 0.001 * obj.lambdaMin;
    
    iteration = 2;
    initialConePositions  = squeeze(obj.latticeAdjustmentSteps(iteration, :, :))*1e6;
    conePositions = smoothGridLocalFunction(hFig, params, grid, initialConePositions);
       
end


function conePositions = smoothGridLocalFunction(hFig, params, gridParams,conePositions)

    positionalDiffTolerance = params.latticeAdjustmentPositionalToleranceF ...
        * gridParams.lambdaMin;
    deps = sqrt(eps) * gridParams.lambdaMin;

    deltaT = 0.2;
    dTolerance = params.latticeAdjustmentDelaunayToleranceF * ...
        gridParams.lambdaMin;

    % Initialize convergence
    oldConePositions = inf;
    forceMagnitudes = [];

    % Turn off Delaunay triangularization warning
    warning('off', 'MATLAB:qhullmx:InternalWarning');

    % Number of cones to be adjusted
    conesNum = size(conePositions, 1);

    % Iteratively adjust the cone positions until the forces between nodes
    % (conePositions) reach equlibrium.
    notConverged = true;
    iteration = 0;
    minDistance = [];
    maxDistance = [];
    meanDistance = [];
    minDistanceForEachCone = [];
    maxYDeltaConeDistanceDisplayed = 5.0;
    
    
    target1ConePosition = [5 0];  
    [~, target1ConeIndex] = min(sum((bsxfun(@minus, conePositions, target1ConePosition)).^2, 2));
    target1ConePosition = conePositions(target1ConeIndex,:);
    %target1ConeIndex = nan;
    
    target2ConePosition = [35 0];  
    [~, target2ConeIndex] = min(sum((bsxfun(@minus, conePositions, target2ConePosition)).^2, 2));
    target2ConePosition = conePositions(target2ConeIndex,:);
    
    target3ConePosition = [50 0];  
    [~, target3ConeIndex] = min(sum((bsxfun(@minus, conePositions, target3ConePosition)).^2, 2));
    target3ConePosition = conePositions(target3ConeIndex,:);
    
    
    visualizedConeXrange = [min([target1ConePosition(1) target3ConePosition(1)])-2 max([target1ConePosition(1) target3ConePosition(1)])+2];
    visualizedConeYrange = 0.5*(target1ConePosition(2)+target3ConePosition(2)) + maxYDeltaConeDistanceDisplayed*[-1 1];


    videoOBJ = VideoWriter('MosaicGeneration.mp4', 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 30;
    videoOBJ.Quality = 100;
    videoOBJ.open();
    
    hFig = figure(10); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1000 270]);
    
    ax = subplot('Position', [0.05 0.1 0.93 0.90]);
    set(ax, 'XLim', [-1 55], 'YLim', [-5 5], 'Color', [1 1 1], 'FontSize', 18);
    axis 'equal'
    set(ax, 'XLim', [1 55], 'YLim', [-5 5]);
    ylabel('\it space (microns)', 'FontSize', 24);
    xlabel('\it space (microns)', 'FontSize', 24);
    box on;
    
    while (notConverged) && (iteration <= params.maxGridAdjustmentIterations)
        iteration = iteration + 1;
    
        % compute cone positional diffs
        positionalDiffs = sqrt(sum((conePositions-oldConePositions) .^ 2, 2));
       
        % check if there are any large movements
        if (max(positionalDiffs) > positionalDiffTolerance)
            % save old come positions
            oldConePositions = conePositions;

            % Perform new Delaunay triangulation to determine the updated
            % topology of the truss. To save computing time, we
            % re-triangulate only when we exceed the
            % positionalDiffTolerance
            triangleConeIndices = delaunayn(conePositions);

            % Compute the centroids of all triangles
            centroidPositions = (conePositions(...
                triangleConeIndices(:, 1), :) + conePositions(...
                triangleConeIndices(:, 2), :) + conePositions(...
                triangleConeIndices(:, 3), :)) / 3;

            % Remove centroids outside the desired region by applying the
            % signed distance function
            d = feval(gridParams.domainFunction, centroidPositions, ...
                gridParams.center, gridParams.radius, ...
                gridParams.ellipseAxes);
            triangleConeIndices = triangleConeIndices(d < ...
                gridParams.borderTolerance, :);

            % Create a list of the unique springs (each spring connecting 2
            % cones)
            % triangleConeIndices is an [M x 3] matrix the m-th row
            % contains indices to the 3 cones that define the triangle
            springs = [triangleConeIndices(:, [1, 2]);
                triangleConeIndices(:, [1, 3]);
                triangleConeIndices(:, [2, 3])];
            springs = unique(sort(springs, 2), 'rows');

            % find all springs connected to this cone
            for coneIndex = 1:conesNum
                springIndices{coneIndex} = find(...
                    (springs(:, 1) == coneIndex) | ...
                    (springs(:, 2) == coneIndex));
            end
            
        end % (max(positionalDiffs) > positionalDiffTolerance)
        
        % Compute spring vectors
        springVectors =  conePositions(springs(:, 1), :) - ...
            conePositions(springs(:, 2), :);
        % their centers
        springCenters = (conePositions(springs(:, 1), :) + ...
            conePositions(springs(:, 2), :)) / 2.0;
        % and their lengths
        springLengths = sqrt(sum(springVectors.^2, 2));
        
        
        % Compute desired spring lengths. This is done by evaluating the
        % passed coneDistance function at the spring centers.
        desiredSpringLengthsAbsolute = feval(gridParams.coneSpacingFunction, ...
            springCenters);

        % Normalize spring lengths
        normalizingFactor = sqrt(sum(springLengths .^ 2) / ...
            sum(desiredSpringLengthsAbsolute .^ 2));
        desiredSpringLengths = desiredSpringLengthsAbsolute * normalizingFactor;

        % Compute spring forces, Force(springLengths, desiredSpringLengths)
        % Force(springLengths, desiredSpringLengths) should be positive
        % when springLengths is near the desiredSpringLengths, which can be
        % achieved by choosing desiredSpringLengths slightly larger than
        % the length we actually desire. Here, we set this to be 1.2
        gain = 1.2;
        springForces = max(gain * desiredSpringLengths - springLengths, 0);

        % compute x, y-components of forces on each of the springs
        springForceXYcomponentVectors = springForces ./ springLengths * ...
            [1, 1] .* springVectors;
        
        springForceXYcomponents = abs(springForceXYcomponentVectors);
        
        springForceTimesDisplacementXYcomponents = springForceXYcomponents;
        
        % Compute net forces on each cone
        netForceVectors = zeros(conesNum, 2);
        for coneIndex = 1:conesNum
            % compute net force from all connected springs
            deltaPos = -bsxfun(@minus, springCenters(...
                springIndices{coneIndex}, :), conePositions(coneIndex, :));

            springForceTimesDisplacementXYcomponents(springIndices{coneIndex}, :) = ...
                sign(deltaPos) .* springForceXYcomponents(springIndices{coneIndex}, :);
            
            netForceVectors(coneIndex, :) = sum(sign(deltaPos) .* ...
                springForceXYcomponents(springIndices{coneIndex}, :), 1);
        end
        
        % force at all fixed cone positions must be 0
        % netForceVectors(1:size(fixedConesPositions, 1), :) = 0;
        
        % Save force magnitudes
        % forceMagnitudes(iteration, :) = ...
        %    sqrt(sum(netForceVectors .^ 2, 2)) / gridParams.lambdaMin;
        
        
        
        renderFrame(ax, iteration, conePositions, visualizedConeXrange, visualizedConeYrange, ...
            target1ConeIndex, target2ConeIndex, target3ConeIndex,...
            springIndices, springForceXYcomponentVectors, desiredSpringLengthsAbsolute, ...
            springForceTimesDisplacementXYcomponents, springCenters, springLengths, netForceVectors);
        
        videoOBJ.writeVideo(getframe(hFig));
        
        
        
        % update cone positions according to netForceVectors
        conePositions = conePositions + deltaT * netForceVectors;
        
        % Find any points that lie outside the domain boundary
        d = feval(gridParams.domainFunction, conePositions, ...
            gridParams.center, gridParams.radius, gridParams.ellipseAxes);
        outsideBoundaryIndices = d > 0;
        
        % And project them back to the domain
        if (~isempty(outsideBoundaryIndices))
            % Compute numerical gradient along x-positions
            dXgradient = (feval(gridParams.domainFunction, ...
                [conePositions(outsideBoundaryIndices, 1) + deps, ...
                conePositions(outsideBoundaryIndices, 2)], ...
                gridParams.center, gridParams.radius, ...
                gridParams.ellipseAxes) - d(outsideBoundaryIndices)) / ...
                deps;
            dYgradient = (feval(gridParams.domainFunction, ...
                [conePositions(outsideBoundaryIndices, 1), ...
                conePositions(outsideBoundaryIndices, 2)+deps], ...
                gridParams.center, gridParams.radius, ...
                gridParams.ellipseAxes) - d(outsideBoundaryIndices)) / ...
                deps;

            % Project these points back to boundary
            conePositions(outsideBoundaryIndices, :) = ...
                conePositions(outsideBoundaryIndices, :) - ...
                [d(outsideBoundaryIndices) .* dXgradient, ...
                d(outsideBoundaryIndices) .* dYgradient];
        end

        % Check if all interior nodes move less than dTolerance
        movementAmplitudes = sqrt(sum(deltaT * netForceVectors(...
            d < -gridParams.borderTolerance, :) .^2 , 2));
        if max(movementAmplitudes) < dTolerance, notConverged = false; end

    end % while (notConverged)
    
end

function renderFrame(ax, iteration, conePositions, visualizedConeXrange, visualizedConeYrange, ...
            target1ConeIndex, target2ConeIndex,target3ConeIndex, ...
            springIndices, springForceXYcomponentVectors, desiredSpringLengthsAbsolute, ...
            springForceTimesDisplacementXYcomponents, springCenters, springLengths, netForceVectors)
        
        % Plot the forces around within the visualized region
        % Find the cone indices to be visualized
        targetConeIndices = find(...
            (conePositions(:,1) >= visualizedConeXrange(1)) & ...
            (conePositions(:,1) <= visualizedConeXrange(2)) & ...
            (conePositions(:,2) >= visualizedConeYrange(1)) & ...
            (conePositions(:,2) <= visualizedConeYrange(2)) );
        
        
        % Plot the springs for the targeted cones
        axes(ax);
        cla(ax);
        hold on;
        for kkk = 1:numel(targetConeIndices)
            targetConeIndex = targetConeIndices(kkk);
            targetConePosition = conePositions(targetConeIndex,:);
            targetConeSpringIndices = springIndices{targetConeIndex};
            plotSprings(targetConePosition, targetConeSpringIndices, ...
                springForceXYcomponentVectors, springCenters, springLengths);
        end
        
        % Plot the cones
        plot(conePositions(:,1), conePositions(:,2), 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 26);
       
        % Outline the targered cones
        if (~isnan(target1ConeIndex))
            plot(conePositions(target1ConeIndex,1), conePositions(target1ConeIndex,2), 'ko', 'MarkerFaceColor', [0.8 0.8 0.3], 'MarkerSize', 26);
        end
        if (~isnan(target2ConeIndex))
            plot(conePositions(target2ConeIndex,1), conePositions(target2ConeIndex,2), 'ko', 'MarkerFaceColor', [0.8 0.8 0.3], 'MarkerSize', 26);
        end
        if (~isnan(target3ConeIndex))
            plot(conePositions(target3ConeIndex,1), conePositions(target3ConeIndex,2), 'ko', 'MarkerFaceColor', [0.8 0.8 0.3], 'MarkerSize', 26);
        end
        
        % Plot the forces on the targeted cones
        if (~isnan(target1ConeIndex))
            target1ConePosition = conePositions(target1ConeIndex,:);
            target1ConeSpringIndices = springIndices{target1ConeIndex};
            plotSpringForces(target1ConePosition, target1ConeSpringIndices, ...
                springForceTimesDisplacementXYcomponents, springCenters, springLengths, desiredSpringLengthsAbsolute);
        end
        
        if (~isnan(target2ConeIndex))
            target2ConePosition = conePositions(target2ConeIndex,:);
            target2ConeSpringIndices = springIndices{target2ConeIndex};
            plotSpringForces(target2ConePosition, target2ConeSpringIndices, ...
                springForceTimesDisplacementXYcomponents, springCenters, springLengths, desiredSpringLengthsAbsolute);
        end
        
        if (~isnan(target3ConeIndex))
            target3ConePosition = conePositions(target3ConeIndex,:);
            target3ConeSpringIndices = springIndices{target3ConeIndex};
            plotSpringForces(target3ConePosition, target3ConeSpringIndices, ...
                springForceTimesDisplacementXYcomponents, springCenters, springLengths, desiredSpringLengthsAbsolute);
        end
        
        
        % Plot the net force vectors on the targeted cones
        if (~isnan(target1ConeIndex))
            plotNetForceVectors(target1ConePosition, netForceVectors(target1ConeIndex, :));
        end
        if (~isnan(target2ConeIndex))
            plotNetForceVectors(target2ConePosition, netForceVectors(target2ConeIndex, :));
        end
        if (~isnan(target3ConeIndex))
            plotNetForceVectors(target3ConePosition, netForceVectors(target3ConeIndex, :));
        end
        title(sprintf('iteration: %2.0f', iteration));
        hold off;
        drawnow;
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
    
    hFig = figure(3);
    clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 400 750 300]);

    rowsNum = figureLayout(1);
    colsNum = figureLayout(2);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rowsNum , ...
           'colsNum', colsNum , ...
           'heightMargin', 0.01, ...
           'widthMargin', 0.00, ...
           'leftMargin', 0.03, ...
           'rightMargin', 0.02, ...
           'bottomMargin', 0.04, ...
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



%%%%% FUNCTIONS FOR SMOOTH GRID

function ellipseAxes = determineEllipseAxesLength(mosaicHalfFOVmicrons)
    largestXYspacing = coneSizeReadData('eccentricity', ...
        1e-6 * mosaicHalfFOVmicrons * [1 1], 'angle', [0 90]);
    if (largestXYspacing(1) < largestXYspacing(2))
        xa = 1;
        ya = largestXYspacing(2) / largestXYspacing(1);
    else
        xa = largestXYspacing(1) / largestXYspacing(2);
        ya = 1;
    end
    ellipseAxes = [xa ya] .^ 2;
end
    
function plotSprings(targetConePosition, targetConeSpringIndices, springForceXYcomponentVectors, springCenters, springLengths)
    for connectedSpringIndex = 1:numel(targetConeSpringIndices)
        springIndex = targetConeSpringIndices(connectedSpringIndex);
        directionVector = springForceXYcomponentVectors(springIndex,:) / norm(springForceXYcomponentVectors(springIndex,:)); 
        % Plot the actual spring lengths
        x1 = springCenters(springIndex,1)-0.5*directionVector(1)*springLengths(springIndex);
        x2 = springCenters(springIndex,1)+0.5*directionVector(1)*springLengths(springIndex);
        y1 = springCenters(springIndex,2)-0.5*directionVector(2)*springLengths(springIndex);
        y2 = springCenters(springIndex,2)+0.5*directionVector(2)*springLengths(springIndex);
        if (sqrt(sum((targetConePosition-[x1 y1]).^2))) < (sqrt(sum((targetConePosition-[x2 y2]).^2)))
            dd = targetConePosition-[x1 y1];
        else
            dd = targetConePosition-[x2 y2];
        end
        plot([x1 x2]+dd(1), [y1 y2]+dd(2), '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 2.0);
    end
end


function plotSpringForces(targetConePosition, targetConeSpringIndices, springForceXYcomponentVectors, springCenters, springLengths, desiredSpringLengthsAbsolute)
    for connectedSpringIndex = 1:numel(targetConeSpringIndices)
        springIndex = targetConeSpringIndices(connectedSpringIndex);
        directionVector = springForceXYcomponentVectors(springIndex,:) / norm(springForceXYcomponentVectors(springIndex,:)); 
%         % Plot the actual spring lengths
%         x1 = springCenters(springIndex,1)-0.5*directionVector(1)*springLengths(springIndex);
%         x2 = springCenters(springIndex,1)+0.5*directionVector(1)*springLengths(springIndex);
%         y1 = springCenters(springIndex,2)-0.5*directionVector(2)*springLengths(springIndex);
%         y2 = springCenters(springIndex,2)+0.5*directionVector(2)*springLengths(springIndex);
%         if (sqrt(sum((targetConePosition-[x1 y1]).^2))) < (sqrt(sum((targetConePosition-[x2 y2]).^2)))
%             dd = targetConePosition-[x1 y1];
%         else
%             dd = targetConePosition-[x2 y2];
%         end
%         plot([x1 x2]+dd(1), [y1 y2]+dd(2), 'w-', 'LineWidth', 2.0);

        % Plot the desired sping lengths
        x1 = springCenters(springIndex,1)-0.5*directionVector(1)*desiredSpringLengthsAbsolute(springIndex);
        x2 = springCenters(springIndex,1)+0.5*directionVector(1)*desiredSpringLengthsAbsolute(springIndex);
        y1 = springCenters(springIndex,2)-0.5*directionVector(2)*desiredSpringLengthsAbsolute(springIndex);
        y2 = springCenters(springIndex,2)+0.5*directionVector(2)*desiredSpringLengthsAbsolute(springIndex);
        if (sqrt(sum((targetConePosition-[x1 y1]).^2))) < (sqrt(sum((targetConePosition-[x2 y2]).^2)))
            dd = targetConePosition-[x1 y1];
        else
            dd = targetConePosition-[x2 y2];
        end
        plot([x1 x2]+dd(1), [y1 y2]+dd(2), 'co-', 'MarkerSize', 14, 'MarkerFaceColor', [0.1 0.8 1], 'LineWidth', 1.5);

        % Plot the repulsive force vectors
        plot(...
            springCenters(springIndex,1)+0.35*[-springForceXYcomponentVectors(springIndex,1) springForceXYcomponentVectors(springIndex,1)], ...
            springCenters(springIndex,2)+0.35*[-springForceXYcomponentVectors(springIndex,2) springForceXYcomponentVectors(springIndex,2)], ...
            'r-', 'LineWidth', 6.0);
    end
end

function plotNetForceVectors(targetConePosition, netForceVector)
    % Plot the net force vector
    pos1 = targetConePosition;
    minDesiredNorm = 0.7;
    maxDesiredNorm = 3.0;
    x = norm(netForceVector);
    desiredNorm = minDesiredNorm + (maxDesiredNorm-minDesiredNorm)*(x-0.05)/(1.5-0.05);

    pos2 = targetConePosition + netForceVector/x*desiredNorm;
    plot([pos1(1) pos2(1)], [pos1(2) pos2(2)], 'k-', 'LineWidth', 4.0);
end

function [coneSpacingInMicrons, eccentricitiesInMicrons] = coneSpacingFunction(conePositions)
    eccentricitiesInMicrons = sqrt(sum(conePositions .^ 2, 2));
    eccentricitiesInMeters = eccentricitiesInMicrons * 1e-6;
    angles = atan2(conePositions(:, 2), conePositions(:, 1)) / pi * 180;
    [coneSpacingInMeters, aperture, density] = ...
        coneSizeReadData('eccentricity', eccentricitiesInMeters, ...
        'angle', angles);
    coneSpacingInMicrons = coneSpacingInMeters' * 1e6;
end

function distances = ellipticalDomainFunction(conePositions, center, radius, ellipseAxes)

    %  points with positive distance will be excluded
    xx = conePositions(:, 1) - center(1);
    yy = conePositions(:, 2) - center(2);
    radii = sqrt((xx / ellipseAxes(1)) .^ 2 + (yy / ellipseAxes(2)) .^ 2);
    distances = -(radius - radii);
end

