function generateFig1Components
% Generate the components for Figure 1 from Cottaris et al 2019(A)
%
% Syntax:
%   generateFig1Components
%
% Description:
%    Generate the necessary components for Figure 1 from the Cottaris et al
%    2019(A) paper.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    [rootPath, ~] = fileparts(which(mfilename));
    rootPath = strrep(rootPath, 'script', 'isetbio_resources');

    generateDisplayFig(rootPath);

    mosaicFOV = 0.6;
    theConeMosaic = generateMosaicFig(rootPath, mosaicFOV);

    sceneFOV = 20;
    visualizedSceneFraction = [0.6 * 1.1 0.6];
    theScene = generateSceneFig(rootPath, sceneFOV, ...
        visualizedSceneFraction, mosaicFOV);
    pause

    opticsModel = 'ThibosDefaultSubject3MMPupil';
    pupilDiamMm = 3.0;
    theOI = generateOpticsImageFig(opticsModel, pupilDiamMm, sceneFOV);

    theOI = oiCompute(theOI, theScene);
    theIsomerizations = theConeMosaic.compute(theOI, 'currentFlag', false);

    generateOpticalImageFig(theOI, ...
        sceneFOV * visualizedSceneFraction, mosaicFOV);
    generateIsomerizationsImageFig(theConeMosaic, theIsomerizations);

end

function generateDisplayFig(rootPath)
% Generate the display figure
%
% Syntax:
%   generateDisplayFig(rootPath)
%
% Description:
%    Generate the required display figure.
%
% Inputs:
%    rootPath - String. The root path to the display information.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    load(fullfile(rootPath, 'display.mat'), 'display');
    hFig = figure(4);
    clf
    set(hFig, 'Position', [10 10 426 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.15 0.96 0.85]);

    wave = [display.wave(1)];
    red = [0];
    green = [0];
    blue = [0];
    for k = 1:numel(display.wave)
        wave = [wave display.wave(k)];
        red = [red display.spd(k, 1)];
        green = [green display.spd(k, 2)];
        blue = [blue display.spd(k, 3)];
        if (k < numel(display.wave))
            wave = [wave display.wave(k + 1)];
            red = [red display.spd(k, 1)];
            green = [green display.spd(k, 2)];
            blue = [blue display.spd(k, 3)];
        end
    end

    area(ax, wave, red, 'EdgeColor', [1 0.0 0.0], 'LineWidth', 1.5, ...
        'FaceColor', [1 0.5 0.5], 'FaceAlpha', .5, 'EdgeAlpha', .8);
    hold(ax, 'on');
    area(ax, wave, green, 'EdgeColor', [0 1.0 0.0], 'LineWidth', 1.5, ...
        'FaceColor', [0.5 1.0 0.5], 'FaceAlpha', .5, 'EdgeAlpha', .8);
    area(ax, wave, blue, 'EdgeColor', [0 0.0 1.0], 'LineWidth', 1.5, ...
        'FaceColor', [0.5 0.5 1], 'FaceAlpha', .5, 'EdgeAlpha', .8);
    grid on;
    set(gca, 'XLim', [display.wave(1) display.wave(end)]);
    set(gca, 'XTick', [400:50:800], 'YTick', [], 'YColor', 'none');
    xlabel('wavelength, \lambda (nm)');
    set(gca, 'FontSize', 20);
    box off
    pwd
    NicePlot.exportFigToPDF('componentFigs/DisplayComponent.pdf', ... 
        hFig, 300);
end

function generateIsomerizationsImageFig(theConeMosaic, theIsomerizations)
% Generate the isomerizations for the image figure
%
% Syntax:
%   generateIsomerizationsImageFig(theConeMosaic, theIsomerizations)
%
% Description:
%    Generate the required isomerizations for the image figure.
%
% Inputs:
%    theConeMosaic     - Struct. A cone mosaic structure.
%    theIsomerizations - 
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    visualizedActivationPattern = theIsomerizations;
    signalRange = [0 max(visualizedActivationPattern(:))];

    hFig = figure(3);
    clf
    set(hFig, 'Position', [10 10 426 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    activationLUT = gray(1024);
    titleForColorBar = 'R*/5 msec';
    backgroundColor = [0 0 0];
    theConeMosaic.integrationTime
    theConeMosaic.renderActivationMap(ax, visualizedActivationPattern, ...
        'signalRange', signalRange, ...
        'visualizedConeAperture', 'geometricArea', ...
        'mapType', 'modulated disks', ...
        'showColorBar', true, ...
        'labelColorBarTicks', true, ...
        'titleForColorBar', titleForColorBar, ...
        'colorMap', activationLUT, ...
        'backgroundColor', backgroundColor);
     ylabel(ax, '');
     set(ax, 'XTickLabels', {});
     set(ax, 'YTickLabels', {});
     NicePlot.exportFigToPDF(...
         'componentFigs/MeanIsomerizationsComponent.pdf', hFig, 300);
end

function generateOpticalImageFig(theOI, xyRange, mosaicFOV)
% Generate the optical image figure
%
% Syntax:
%   generateOpticalImageFig(theOI, xyRange, mosaicFOV)
%
% Description:
%    Generate the figure containing the optical image.
%
% Inputs:
%    theOI     - Struct. An optical image structure.
%    xyRange   - Vector. A 1x2 vector containing the X and Y ranges.
%    mosaicFOV - Numeric. The mosaic's field of view, in degrees.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    rgbImage = oiGet(theOI, 'rgb image');
    sampleSizeDegs = oiGet(theOI, 'wangular resolution');
    xSupport = 1:size(rgbImage, 2);
    xSupport = xSupport * sampleSizeDegs;
    xSupport = xSupport - mean(xSupport);
    ySupport = 1:size(rgbImage, 2);
    ySupport = ySupport * sampleSizeDegs;
    ySupport = ySupport - mean(ySupport);

    hFig = figure(2);
    clf;
    set(hFig, 'Position', [10 10 468 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    image(ax, xSupport , ySupport, rgbImage);
    hold(ax, 'on');
    x = 0.5 * mosaicFOV * [-1 1 1 -1 -1];
    y = 0.5 * mosaicFOV * [-1 -1 1 1 -1];
    plot(ax, x, y, 'k-', 'LineWidth', 1.0);
    hold(ax, 'off');
    axis(ax, 'xy');
    axis(ax, 'image');

    tickInterval = 0.2;
    set(ax, 'XLim', xyRange(1) / 2 * [-1 1], ...
        'YLim', xyRange(2) / 2 * [-1 1], ...
         'XTick', -15:tickInterval:15, 'YTick', -15:tickInterval:15, ...
         'YTickLabels', {}, 'YTickLabels', {});
    grid on;
    box on;
    set(ax, 'FontSize', 12);
    %xlabel('space (degs)');
    %ylabel('space (degs)');
    NicePlot.exportFigToPNG('componentFigs/OpticalImageComponent.png', ...
        hFig, 300);

    wave = oiGet(theOI, 'wave');
    photons = oiGet(theOI, 'photons');
    maxPhotons = max(photons(:));
    minPhotons = min(photons(:));
    visualizedWaves = [450 480 510 530 550 580 610 640 670];
    for k = 1:numel(visualizedWaves)
        [~, idx] = min(abs(wave-visualizedWaves(k)));
        thePhotons = squeeze(photons(:, :, idx));
        hFig = figure(2999 + k);
        clf;
        set(hFig, 'Position', [10 10 468 420], 'Color', [1 1 1]);
        ax = subplot('Position', [0.02 0.02 0.96 0.96]);
        imagesc(ax, xSupport, ySupport, ...
            (thePhotons - minPhotons) / (maxPhotons - minPhotons), [0 1]);
        hold(ax, 'on');
        x = 0.5 * mosaicFOV * [-1 1 1 -1 -1];
        y = 0.5 * mosaicFOV * [-1 -1 1 1 -1];
        plot(ax, x, y, 'k-', 'LineWidth', 1.0);
        hold(ax, 'off');
        axis(ax, 'image');
        axis(ax, 'xy');
        set(ax, 'XLim', xyRange(1) / 2 * [-1 1], ...
            'YLim', xyRange(2) / 2 * [-1 1], ...
         'XTick', -15:tickInterval:15, 'YTick', -15:tickInterval:15, ...
         'XTickLabels', {}, 'YTickLabels', {});
        colormap(gray(1024));
        NicePlot.exportFigToPNG(sprintf(...
            'componentFigs/PhotonsOIComponent%d.png', wave(idx)), ...
            hFig, 300);
    end

    micronsPerDegree = 300;
    visualizePSFfromOI(theOI, micronsPerDegree, ...
        'colormapToUse', gray(1024), ...
        'visualizedWavelengths', [450 480 510 530 550 580 610 640 670], ...
        'rows', 3, 'cols', 3, ...
        'labelLastPSF', false, ...
        'displayWavelengthInTitle', ~false);

end

function theScene = generateSceneFig(rootPath, sceneFOV, ...
    visualizedSceneFraction, mosaicFOV)
% Generate the scene figure
%
% Syntax:
%   theScene = generateSceneFig(rootPath, sceneFOV, ...
%       visualizedSceneFraction, mosaicFOV)
%
% Description:
%    Generate the scene figure.
%
% Inputs:
%    rootPath  - String. The fullfile path to the scene information.
%    sceneFOV  - Numeric. The scene's field of view, in degrees.
%    visualizedSceneFraction
%              - Numeric. The fraction of the scene actually being shown.
%    mosaicFOV - Numeric. The mosaic's field of view, in degrees.
%
% Outputs:
%    theScene  - Struct. A scene structure.
%
% Optional key/value pairs:
%    None.
%

    load(fullfile(rootPath, 'scenes.mat'));
    theScene = scenesList{1};

    theScene = sceneSet(theScene, 'h fov', sceneFOV);
    sceneGet(theScene, 'wAngular')

    save(fullfile('componentFigs', 'scene.mat'), 'theScene');

    hFig = figure(1);
    clf
    set(hFig, 'Position', [10 10 468 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    sceneSize = sceneGet(theScene, 'size');
    rows = sceneSize(1);
    cols = sceneSize(2);
    sceneGet(theScene, 'width', 'm')
    xSupport = (1:cols);
    xSupport = xSupport - mean(xSupport);
    xSupport = 0.5 * xSupport / max(abs(xSupport));
    xSupport = xSupport * sceneGet(theScene, 'wAngular');
    ySupport = (1:rows);
    ySupport = ySupport - mean(ySupport);
    ySupport = 0.5 * ySupport / max(abs(ySupport));
    ySupport = ySupport * sceneGet(theScene, 'hAngular');

    imagesc(ax, xSupport, ySupport, sceneGet(theScene, 'rgb image'));
    hold(ax, 'on');
    x = 0.5 * mosaicFOV * [-1 1 1 -1 -1];
    y = 0.5 * mosaicFOV * [-1 -1 1 1 -1];
    plot(ax, x, y, 'k-', 'LineWidth', 1.0);
    hold(ax, 'off');

    axis(ax, 'image');
    axis(ax, 'xy');
    tickInterval = 0.2;
    xyRange = sceneFOV * visualizedSceneFraction;
    set(ax, 'XLim', xyRange(1) / 2 * [-1 1], ...
        'YLim', xyRange(2) / 2 * [-1 1], ...
        'XTick', -15:tickInterval:15, 'YTick', -15:tickInterval:15, ...
        'XTickLabels', {}, 'YTickLabels', {});
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontSize', 12);
    % xlabel('space (degs)');
    % ylabel('space (degs)');
    NicePlot.exportFigToPNG('componentFigs/SceneComponent.png', hFig, 300);

    wave = sceneGet(theScene, 'wave');
    photons = sceneGet(theScene, 'photons');
    maxPhotons = max(photons(:));
    minPhotons = min(photons(:));
    visualizedWaves = [450 480 510 530 550 580 610 640 670];
    for k = 1:numel(visualizedWaves)
        [~, idx] = min(abs(wave - visualizedWaves(k)));
        thePhotons = squeeze(photons(:, :, idx));
        hFig = figure(1999 + k);
        clf;
        set(hFig, 'Position', [10 10 468 420], 'Color', [1 1 1]);
        ax = subplot('Position', [0.02 0.02 0.96 0.96]);
        imagesc(ax, xSupport, ySupport, ...
            (thePhotons - minPhotons) / (maxPhotons - minPhotons), [0 1]);
        hold(ax, 'on');
        x = 0.5*mosaicFOV * [-1 1 1 -1 -1];
        y = 0.5*mosaicFOV * [-1 -1 1 1 -1];
        plot(ax, x, y, 'k-', 'LineWidth', 1.0);
        hold(ax, 'off');
        axis(ax, 'image');
        axis(ax, 'xy');
        set(ax, 'XLim', xyRange(1) / 2 * [-1 1], ...
            'YLim', xyRange(2) / 2 * [-1 1], ...
            'XTick', -15:tickInterval:15, 'YTick', -15:tickInterval:15, ...
            'XTickLabels', {}, 'YTickLabels', {});
        colormap(gray(1024));
        NicePlot.exportFigToPNG(sprintf(...
            'componentFigs/PhotonsSceneComponent%d.png', wave(idx)), ...
            hFig, 300);
    end
end

function theConeMosaic = generateMosaicFig(rootPath, mosaicFOV)
% Generate the mosaic figure
%
% Syntax:
%   theConeMosaic = generateMosaicFig(rootPath, mosaicFOV)
%
% Description:
%    Generate the figure for the cone mosaic.
%
% Inputs:
%    rootPath      - String. The fullfile path to the cone mosaic data.
%    mosaicFOV     - Numeric. The mosaic field of view, in degrees.
%
% Outputs:
%    theConeMosaic - Struct. A cone mosaic structure.
%
% Optional key/value pairs:
%    None.
%

    load(fullfile(rootPath, sprintf(...
        'coneMosaic_%1.2fdegFOV.mat', mosaicFOV)), 'theConeMosaic');

    hFig = figure(4);
    clf
    set(hFig, 'Position', [10 10 426 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);

    theConeMosaic.visualizeGrid('axesHandle', ax, ...
        'visualizedConeAperture', 'geometricArea');
    set(ax, 'XColor', 'none', 'YColor', 'none');
    set(ax, 'XTickLabels', {});
    set(ax, 'YTickLabels', {});
    NicePlot.exportFigToPDF('componentFigs/MosaicComponent.pdf', ...
        hFig, 300);
end

function theOI = ...
    generateOpticsImageFig(opticsModel, pupilDiamMm, horizontalFOV)
% Generate an image figure for the optics model structure.
%
% Syntax:
%   theOI = generateOpticsImageFig(opticsModel, pupilDiamMm, horizontalFOV)
%
% Description:
%    Generate an image figure for the Optics model.
%
% Inputs:
%    opticsModel   - Struct. A structure for the optics model.
%    pupilDiamMm   - Numeric. The pupil diameter in millimeters.
%    horizontalFOV - Numeric. The horizontal field of view in degrees.
%
% Outputs:
%    theOI         - Struct. The optical image structure.
%
% Optional key/value pairs:
%    None.
%

    oiParams = struct('opticsModel', opticsModel, ...
        'wavefrontSpatialSamples', 301, 'pupilDiamMm', pupilDiamMm, ...
        'umPerDegree', 300);
    theOI = oiWithCustomOptics(oiParams.opticsModel, ...
        oiParams.wavefrontSpatialSamples, oiParams.pupilDiamMm, ...
        oiParams.umPerDegree);

    % Set the FOV
    theOI = oiSet(theOI, 'h fov', horizontalFOV);

    % Set the fNumber
    focalLength = oiGet(theOI, 'distance');
    desiredFNumber = focalLength / (oiParams.pupilDiamMm / 1000);
    theOI  = oiSet(theOI , 'optics fnumber', desiredFNumber);
end
