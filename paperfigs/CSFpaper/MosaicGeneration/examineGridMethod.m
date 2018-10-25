function examineGridMethod
    load('s.mat')
    params.latticeAdjustmentPositionalToleranceF = s.positionalToleranceF/8;
    params.latticeAdjustmentDelaunayToleranceF = s.DelaunayToleranceF/8;
    params.maxGridAdjustmentIterations = s.maxGridAdjustmentIterations*8;
    params.saveLatticeAdjustmentProgression = true;
    params.latticeAdjustmentSteps = [];
    
    conePositions = s.conePositions;
    gridParams = s.gridParams;

    
    size(conePositions)
    hFig = figure(1);
    set(hFig, 'Position', [10 10 1000 1000], 'Color', [1 1 1]);
    
    
    conePositions = smoothGrid(hFig, params, gridParams, conePositions);
end

function [minDistance, distanceToClosestNeighbor] = findMinDistance(conePositions)
    D = squareform(pdist(conePositions));
    % Make diagonal points (which have zero distance - distance from a
    % point to itself) equal to nan
    D(logical(eye(size(D)))) = nan;
    minDistance = min(D(:));
    conesNum = size(D,1);
    distanceToClosestNeighbor = zeros(1,conesNum);
    for coneIndex = 1:conesNum
        distancesToAllCones = squeeze(D(coneIndex,:));
        distanceToClosestNeighbor(1,coneIndex) = min(distancesToAllCones);
    end
    
    %[sortedD, indices] = sort(D(:));
    %[min(sortedD(:)) max(sortedD(:))]
end


function conePositions = smoothGrid(hFig, obj, gridParams,conePositions)


    positionalDiffTolerance = obj.latticeAdjustmentPositionalToleranceF ...
        * gridParams.lambdaMin;
    deps = sqrt(eps) * gridParams.lambdaMin;

    deltaT = 0.2;
    dTolerance = obj.latticeAdjustmentDelaunayToleranceF * ...
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
    tic
    minDistance = [];
    maxDistance = [];
    meanDistance = [];
    minDistanceForEachCone = [];
    maxDeltaConeDistanceDisplayed = 6;
    
    target1ConePosition = [22 30];  
    [~, target1ConeIndex] = min(sum((bsxfun(@minus, conePositions, target1ConePosition)).^2, 2));
    target1ConePosition = conePositions(target1ConeIndex,:);
    %target1ConeIndex = nan;
    
    target2ConePosition = [24 33];  
    [~, target2ConeIndex] = min(sum((bsxfun(@minus, conePositions, target2ConePosition)).^2, 2));
    target2ConePosition = conePositions(target2ConeIndex,:);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 2, ...
                'colsNum', 2, ...
                'heightMargin', 0.07, ...
                'widthMargin', 0.06, ...
                'leftMargin', 0.06, ...
                'rightMargin', 0.001, ...
                'bottomMargin', 0.05, ...
                'topMargin', 0.001);
            
    videoOBJ = VideoWriter('MosaicGeneration.mp4', 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 30;
    videoOBJ.Quality = 100;
    videoOBJ.open();
                
    while (notConverged) && (iteration <= obj.maxGridAdjustmentIterations)
        iteration = iteration + 1;
        [~, distanceToClosestNeighbor(iteration,:)] = findMinDistance(conePositions);
        
        figure(hFig); clf;
        subplot('Position', subplotPosVectors(1,1).v);
        plot(conePositions(:,1), conePositions(:,2), 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 9);
        hold on;
        if (~isnan(target1ConeIndex))
            plot(conePositions(target1ConeIndex,1), conePositions(target1ConeIndex,2), 'ro', 'MarkerFaceColor', [0.9 0.5 0.5], 'MarkerSize', 9);
        end
        if (~isnan(target2ConeIndex))
            plot(conePositions(target2ConeIndex,1), conePositions(target2ConeIndex,2), 'go', 'MarkerFaceColor', [0.5 0.9 0.5], 'MarkerSize', 9);
        end
        set(gca, 'XLim', 40*[-1 1], 'YLim', 40*[-1 1], 'Color', [0 0 0], 'FontSize', 18);
        axis 'square'
        xlabel('space (microns)');
        ylabel('space (microns)');
        
        subplot('Position', subplotPosVectors(1,2).v);
        visualizedConeXrange = 0.5*(target1ConePosition(1)+target2ConePosition(1)) + maxDeltaConeDistanceDisplayed*[-1 1];
        visualizedConeYrange = 0.5*(target1ConePosition(2)+target2ConePosition(2)) + maxDeltaConeDistanceDisplayed*[-1 1];
        set(gca, 'XLim', visualizedConeXrange, ...
                 'YLim', visualizedConeYrange, 'Color', [0 0 0], 'FontSize', 18);
        axis 'square'
        xlabel('space (microns)');
        
        subplot('Position', subplotPosVectors(2,1).v);
        plot(1:size(distanceToClosestNeighbor,2), squeeze(distanceToClosestNeighbor(iteration,:)), 'ko', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 5);
        set(gca, 'XLim', [1 size(distanceToClosestNeighbor,2)], 'YLim', [0 5], 'YTick', [0:5]);
        set(gca, 'FontSize', 18);
        xlabel('cone index');
        ylabel('distance to closest neighbor (microns)');
        
        grid on;
        axis 'square'
        
        subplot('Position', subplotPosVectors(2,2).v);
        minDistance(iteration) = min(squeeze(distanceToClosestNeighbor(iteration,:)));
        maxDistance(iteration) = max(squeeze(distanceToClosestNeighbor(iteration,:)));
        meanDistance(iteration) = mean(squeeze(distanceToClosestNeighbor(iteration,:)));
        fprintf('[%d]: Min distance b/n all points: %2.2f\n', iteration, minDistance(iteration));
        if (iteration < 30)
            markerSize = 12;
        elseif (iteration < 50)
            markerSize = 11;
        elseif (iteration < 100)
            markerSize = 10;
        elseif (iteration < 200)
            markerSize = 9;
        elseif (iteration < 400)
            markerSize = 8;
        elseif (iteration < 800)
            markerSize = 6;
        else
            markerSize = 4;
        end
        plot(1:iteration, minDistance, 'ro-', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', markerSize, 'LineWidth', 1.5);
        hold on;
        plot(1:iteration, maxDistance, 'bo-',  'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', markerSize, 'LineWidth', 1.5);
        plot(1:iteration, meanDistance, 'ko-',  'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', markerSize, 'LineWidth', 1.5);
        hold off;
        set(gca, 'XTick', [0:20:1000], 'YLim', [0 5], 'YTick', [0:5], 'XLim', [0 max([10 iteration])]);
        if (iteration <= 10)
            set(gca, 'XTick', [0:1:10]);
        elseif (iteration <= 50)
            set(gca, 'XTick', [0:5:50]);
        elseif (iteration <= 100)
            set(gca, 'XTick', [0:10:100]);
        elseif (iteration <= 500)
            set(gca, 'XTick', [0:50:500]);
        elseif (iteration <= 1000)
            set(gca, 'XTick', [0:100:iteration]);
        else
             set(gca, 'XTick', [0:250:iteration]);
        end
        
        set(gca, 'FontSize', 18);
        grid on;
        legend({'min', 'max', 'mean'});
        xlabel('iteration #');
        ylabel('distance to closest neighbor (microns)');
        axis 'square'
        

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
        
        
        subplot('Position', subplotPosVectors(1,2).v);
        hold on;
        
        % Plot the springs
        targetConeIndices = find(...
            (conePositions(:,1) >= visualizedConeXrange(1)) & ...
            (conePositions(:,1) <= visualizedConeXrange(2)) & ...
            (conePositions(:,2) >= visualizedConeYrange(1)) & ...
            (conePositions(:,2) <= visualizedConeYrange(2)) );
        
        for kkk = 1:numel(targetConeIndices)
            targetConeIndex = targetConeIndices(kkk);
            targetConePosition = conePositions(targetConeIndex,:);
            targetConeSpringIndices = springIndices{targetConeIndex};
            plotSprings(targetConePosition, targetConeSpringIndices, ...
                springForceXYcomponentVectors, springCenters, springLengths);
        end

        % Plot the cones
        plot(conePositions(:,1), conePositions(:,2), 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 26);
        if (~isnan(target1ConeIndex))
            plot(conePositions(target1ConeIndex,1), conePositions(target1ConeIndex,2), 'ro', 'MarkerFaceColor', [0.9 0.5 0.5], 'MarkerSize', 26);
        end
        if (~isnan(target2ConeIndex))
            plot(conePositions(target2ConeIndex,1), conePositions(target2ConeIndex,2), 'go', 'MarkerFaceColor', [0.4 0.8 0.4], 'MarkerSize', 26);
        end
        % Plot the forces
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
        
        if (~isnan(target1ConeIndex))
            plotNetForceVectors(target1ConePosition, netForceVectors(target1ConeIndex, :));
        end
        if (~isnan(target2ConeIndex))
            plotNetForceVectors(target2ConePosition, netForceVectors(target2ConeIndex, :));
        end
        
        drawnow;
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

        if (obj.saveLatticeAdjustmentProgression)
            obj.latticeAdjustmentSteps(...
                size(obj.latticeAdjustmentSteps, 1) + 1, :, :) = ...
                conePositions * 1e-6;
        end 
    end % while (notConverged) && (iteration < obj.maxGridAdjustmentIterations)
    
    videoOBJ.close();
    
    fprintf('\nHex grid smoothing finished in %2.1f seconds.', toc);
    if (iteration > obj.maxGridAdjustmentIterations) 
        fprintf('\nDid not converge, but exceeded max number of iterations (%d).', ...
            obj.maxGridAdjustmentIterations);
        fprintf('\nMax(movement) in last iteration: %2.6f, Tolerange: %2.6f\n', ...
            max(movementAmplitudes), obj.latticeAdjustmentDelaunayToleranceF);
    else
        fprintf('Converged after %d iterations.\n', iteration);
    end
    
    % Turn back on Delaunay triangularization warning
    warning('on', 'MATLAB:qhullmx:InternalWarning');
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
        plot([x1 x2]+dd(1), [y1 y2]+dd(2), 'w-', 'LineWidth', 2.0);
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
        plot([x1 x2]+dd(1), [y1 y2]+dd(2), 'co-', 'MarkerSize', 12, 'MarkerFaceColor', [0.1 0.8 .8], 'LineWidth', 2.0);

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
    plot([pos1(1) pos2(1)], [pos1(2) pos2(2)], 'y-', 'LineWidth', 4.0);
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
