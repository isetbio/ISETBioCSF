function exploreMosaic

    [rootPath0,~] = fileparts(which(mfilename));
    rootPath = strrep(rootPath0, 'sideprojects/MosaicGenerator', 'FVM2018Scripts/resources');
    
    mosaicFOV = 20;
    mosaicFileName = sprintf('ConeMosaic_%2.1fDegs.mat', mosaicFOV);
    %load(mosaicFileName, 'theConeMosaic');
    load('/Volumes/SamsungT3/MATLAB/projects/ISETBioLiveScript/toolbox/resources/ConeMosaic20Degs.mat');
    pause
    
    roiRect = struct(...
        'xo', 6, 'yo',.0, ...
        'width', 3.0, 'height', 2.0);
    
    roiRect2 = struct(...
        'xo', -1.5, 'yo', -1.0, ...
        'width', 2.0, 'height', 3.0);
    
    theDisplay = generateTheDisplay();
    theScene = generateTheScene(rootPath,theDisplay, mosaicFOV*1.5);
    theOI = generateTheOpticalImage(theScene);
    theOI = oiCompute(theOI, theScene);
    save('OpticalImage20Degs.mat', 'theOI');
    pause
    
    hFig = figure(2); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1620 840]);
    subplot(2,5,3);
    visualizeOI(theOI, mosaicFOV, roiRect);

    ax = subplot(2,5,[1 2 6 7]);
    theMosaic.visualizeGrid('axesHandle', ax, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'geometricArea', ...
        'ticksInVisualDegs', true, ...
        'tickInc', 5.0, ...
        'backgroundColor', [1 1 1]);
    drawnow;
    
    % Clip mosaic
    theMosaic.clipWithRect(roiRect);
    
    % Compute isomerizations
    theIsomerizations =  theMosaic.compute(theOI, 'currentFlag',false);
    maxIsometizations1 = max(theIsomerizations(:));
    
    ax = subplot(2,5,4);
    theMosaic.visualizeGrid('axesHandle', ax, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'geometricArea', ...
        'visualizedFOV', roiRect, ...
        'tickInc', 1.0, ...
        'ticksInVisualDegs', true, ...
        'backgroundColor', [1 1 1]);
    
    ax = subplot(2,5,5);
    theMosaic.renderActivationMap(ax, theIsomerizations, ...
        'mapType', 'modulated disks', ...
        'visualizedConeAperture', 'geometricArea', ...
        'visualizedFOV', roiRect, ...
        'tickInc', 1.0, ...
        'colorMap', gray(1024), ...
        'backgroundColor', [0 0 0]);
    set(gca, 'XTickLabel', {}, 'YTickLabel', {});
    ylabel('');
    
    % Reload mosaic
    load(mosaicFileName);
    
    % Clip mosaic with rect 2
    theMosaic.clipWithRect(roiRect2);
    
    
   % visualizeScene(theScene, mosaicFOV, roiRect2);
    
    subplot(2,5,8);
    visualizeOI(theOI, mosaicFOV, roiRect2);

    % Compute isomerizations
    theIsomerizations =  theMosaic.compute(theOI, 'currentFlag',false);
    maxIsometizations2 = max(theIsomerizations(:));
    
    ax = subplot(2,5,9);
    theMosaic.visualizeGrid('axesHandle', ax, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'geometricArea', ...
        'visualizedFOV', roiRect2, ...
        'tickInc', 1.0, ...
        'ticksInVisualDegs', true, ...
        'backgroundColor', [1 1 1]);
    
    ax = subplot(2,5,10);
    theMosaic.renderActivationMap(ax, theIsomerizations, ...
        'mapType', 'modulated disks', ...
        'visualizedConeAperture', 'geometricArea', ...
        'visualizedFOV', roiRect2, ...
        'tickInc', 1.0, ...
        'colorMap', gray(1024), ...
        'backgroundColor', [0 0 0]);
    set(gca, 'XTickLabel', {}, 'YTickLabel', {});
    ylabel('');
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
    set(gca, 'XTick', -10:5:10, 'YTick', -10:5:10);
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
    set(gca, 'XTick', -10:5:10, 'YTick', -10:5:10);
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