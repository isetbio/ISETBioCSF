function visualizeLargeMosaic()
    load('HexConeMosaic15DegsBest.mat', 'theMosaic');
    
    % Visualize only part of the mosaic X(degs): [0 6], Y(degs): [-0.25 0.25]
    visualizedFOV = struct('xo', 3.0, 'yo', 0.0, 'width', 6, 'height', 0.5);
    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 2800 600]);
    axesHandle  = subplot('Position', [0.02 0.5 0.97 0.45]);
    theMosaic.visualizeGrid('axesHandle', axesHandle, 'ticksInVisualDegs', true, 'visualizedFOV', visualizedFOV);
end