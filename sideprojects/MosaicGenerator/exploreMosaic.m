function exploreMosaic

    [rootPath0,~] = fileparts(which(mfilename));
    rootPath = strrep(rootPath0, 'sideprojects/MosaicGenerator', 'FVM2018Scripts/resources');
    
    mosaicFOV = 10;
    posTolerance = 1.5;
    mosaicFileName = sprintf('ConeMosaic_%2.1fDegs_PosTolerance%2.2f.mat', mosaicFOV,posTolerance);
    load(mosaicFileName);
    
    
    roiRect = struct(...
        'xo', -0.3, 'yo', 0.3, ...
        'width', 1.0, 'height', 0.6);
    
    roiRect2 = struct(...
        'xo', 1.5, 'yo', -0.3, ...
        'width', 1.0, 'height', 1.0);
    
    theDisplay = generateTheDisplay();
    theScene = generateTheScene(rootPath,theDisplay, mosaicFOV*1.5);
    theOI = generateTheOpticalImage(theScene);
    theOI = oiCompute(theOI, theScene);
    

    hFig = figure(2); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1693 1344]);
    subplot(4,3,1);
    visualizeScene(theScene, mosaicFOV, roiRect);
    
    subplot(4,3,2)
    visualizeOI(theOI, mosaicFOV, roiRect);

    ax = subplot(4,3,3);
    theMosaic.visualizeGrid('axesHandle', ax, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'geometricArea', ...
        'ticksInVisualDegs', true, ...
        'backgroundColor', [1 1 1]);
    drawnow;
    
    % Clip mosaic
    theMosaic.clipWithRect(roiRect);
    
    % Compute isomerizations
    theIsomerizations =  theMosaic.compute(theOI, 'currentFlag',false);
    maxIsometizations1 = max(theIsomerizations(:))
    
    ax = subplot(4,3,4);
    theMosaic.visualizeGrid('axesHandle', ax, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'geometricArea', ...
        'visualizedFOV', roiRect, ...
        'ticksInVisualDegs', true, ...
        'backgroundColor', [1 1 1]);
    
    ax = subplot(4,3,5);
    theMosaic.renderActivationMap(ax, theIsomerizations, ...
        'mapType', 'modulated disks', ...
        'visualizedConeAperture', 'geometricArea', ...
        'visualizedFOV', roiRect, ...
        'colorMap', gray(1024), ...
        'backgroundColor', [0 0 0]);
    
    % Reload mosaic
    load(mosaicFileName);
    
    % Clip mosaic with rect 2
    theMosaic.clipWithRect(roiRect2);
    
    
    subplot(4,3,7);
    visualizeScene(theScene, mosaicFOV, roiRect2);
    
    subplot(4,3,8)
    visualizeOI(theOI, mosaicFOV, roiRect2);

    % Compute isomerizations
    theIsomerizations =  theMosaic.compute(theOI, 'currentFlag',false);
    maxIsometizations2 = max(theIsomerizations(:))
    
    ax = subplot(4,3,10);
    theMosaic.visualizeGrid('axesHandle', ax, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'geometricArea', ...
        'visualizedFOV', roiRect2, ...
        'ticksInVisualDegs', true, ...
        'backgroundColor', [1 1 1]);
    
    ax = subplot(4,3,11);
    theMosaic.renderActivationMap(ax, theIsomerizations, ...
        'mapType', 'modulated disks', ...
        'visualizedConeAperture', 'geometricArea', ...
        'visualizedFOV', roiRect2, ...
        'colorMap', gray(1024), ...
        'backgroundColor', [0 0 0]);
    
end


function visualizeScene(theScene, mosaicFOV, roiRect)
    rgbImage = sceneGet(theScene, 'rgb image');
    sceneSize = sceneGet(theScene, 'size');
    rows = sceneSize(1); cols = sceneSize(2);
    xSupport = (1:cols);
    xSupport = xSupport - mean(xSupport);
    xSupport = 0.5*xSupport / max(abs(xSupport));
    sceneWidth = sceneGet(theScene, 'wAngular');
    sceneHeight = sceneGet(theScene, 'hAngular');
    sceneFOV = [sceneWidth(1) sceneHeight(1)];
    xSupport = xSupport * sceneWidth(1);
    ySupport = (1:rows);
    ySupport = ySupport - mean(ySupport);
    ySupport = 0.5*ySupport / max(abs(ySupport));
    ySupport = ySupport * sceneHeight(1);
    
    roiX = roiRect.xo + roiRect.width/2* [-1 -1 1 1 -1];
    roiY = roiRect.yo + roiRect.height/2*[-1 1 1 -1 -1];
    imagesc(xSupport, ySupport, rgbImage); hold on;
    plot(roiX, roiY, 'r-', 'LineWidth', 1.5);
    hold off
    axis 'xy'
    axis 'image'
    set(gca, 'XLim', 0.5*max(mosaicFOV)*[-1 1], 'YLim', 0.5*max(mosaicFOV)*[-1 1]);
    set(gca, 'FontSize', 18);
    xlabel('space (degs)');
    ylabel('space (degs)');
    
end

function  visualizeOI(theOI, mosaicFOV, roiRect)
    rgbImage = oiGet(theOI, 'rgb image');
    sampleSizeDegs = oiGet(theOI, 'wangular resolution');
    xSupport = 1:size(rgbImage,2);
    xSupport = xSupport * sampleSizeDegs;
    xSupport = xSupport - mean(xSupport);
    ySupport = 1:size(rgbImage,1);
    ySupport = ySupport * sampleSizeDegs;
    ySupport = ySupport - mean(ySupport);
    imagesc(xSupport, ySupport, rgbImage);
    roiX = roiRect.xo + roiRect.width/2* [-1 -1 1 1 -1];
    roiY = roiRect.yo + roiRect.height/2*[-1 1 1 -1 -1];
    imagesc(xSupport, ySupport, rgbImage); hold on;
    plot(roiX, roiY, 'r-', 'LineWidth', 1.5);
    hold off
    axis 'xy'
    axis 'image'
    set(gca, 'XLim', 0.5*max(mosaicFOV)*[-1 1], 'YLim', 0.5*max(mosaicFOV)*[-1 1]);
    set(gca, 'FontSize', 18);
    xlabel('space (degs)');
    ylabel('space (degs)');
    
end

function theOI = generateTheOpticalImage(theScene)
    opticsModel = 'ThibosDefaultSubject3MMPupil'; 
    pupilDiamMm = 3.0;
    
    oiParams = struct(...
                'opticsModel', opticsModel, ...
                'wavefrontSpatialSamples', 301, ...
                'pupilDiamMm', pupilDiamMm, ...
                'umPerDegree', 300);
            
    if (strcmp(opticsModel, 'Geisler'))
        theOI = oiCreate('wvf human', oiParams.pupilDiamMm,[],[], oiParams.umPerDegree);
        theOI = ptb.oiSetPtbOptics(theOI,'opticsModel',oiParams.opticsModel);
    else
        theOI = oiWithCustomOptics(oiParams.opticsModel, oiParams.wavefrontSpatialSamples, oiParams.pupilDiamMm, oiParams.umPerDegree);   
    end
    
    % Set the fNumber
    focalLength = oiGet(theOI,'distance');
    desiredFNumber = focalLength/(oiParams.pupilDiamMm/1000);
    theOI  = oiSet(theOI ,'optics fnumber',desiredFNumber);
end

function [theScene, sensorImageRGB] = generateTheScene(rootPath,theDisplay, theMosaicFOV)
%   
    sceneName = 'Blobbie'; % 'Gabor';
    if strcmp(sceneName,'Blobbie')
        rgbSettingsFile = 'PilotC4M4-RGB.mat';  % matte
        rgbSettingsFile = 'PilotC1M4-RGB.mat';  % shiny
        filename = fullfile(rootPath,rgbSettingsFile);
        load(filename, 'sensorImageRGB');
        sensorImageRGB = flipud(sensorImageRGB);
        
        %% Convert to ISETbio scene
        meanLuminance = [];
        theScene = sceneFromFile(sensorImageRGB, 'rgb', meanLuminance, theDisplay);

        meanLuminance = sceneGet(theScene, 'mean luminance');
    else
        stimParams = struct(...
            'spatialFrequencyCyclesPerDeg', 16, ... % 15 cycles/deg
            'orientationDegs', 0, ...               % 45 degrees
            'phaseDegs', 90, ...                    % spatial phase in degrees
            'sizeDegs', 0.5, ...                    % 0.5 x 0.5 degrees
            'sigmaDegs', 0.2/3, ...                 % sigma of Gaussian envelope
            'contrast', 0.9,...                     % 0.6 Michelson contrast
            'meanLuminanceCdPerM2', 40, ...         % mean luminance
            'pixelsAlongWidthDim', 256, ...         % pixels- width dimension
            'pixelsAlongHeightDim', 256 ...         % pixels- height dimension
            );
        theScene =  generateGaborScene('stimParams', stimParams);
        sensorImageRGB = [];
    end
    
    theScene = sceneAdjustLuminance(theScene, 100);
    meanLuminance = sceneGet(theScene, 'mean luminance');
    % Scale the scene to be a little larger than the mosaic
    sceneSize = sceneGet(theScene, 'size');
    theScene = sceneSet(theScene, 'h fov', theMosaicFOV(1));
end


function theDisplay = generateTheDisplay()
    backgroundParams.backgroundxyY = [0.27 0.30 49.8]';
    backgroundParams.monitorFile = 'CRT-MODEL';
    backgroundParams.leakageLum = 1.0;
    backgroundParams.lumFactor = 1;
    
    theDisplay = displayCreate(backgroundParams.monitorFile);
end