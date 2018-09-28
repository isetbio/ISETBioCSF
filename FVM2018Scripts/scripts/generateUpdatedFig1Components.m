function generateUpdatedFig1Components
    [rootPath0,~] = fileparts(which(mfilename));
    cd(rootPath0);
    rootPath = strrep(rootPath0, 'scripts', 'resources');
    videoOutDir = strrep(rootPath0, 'scripts', 'updatedComponentFigs');
    
    %hFig = runConditionToVisualizePsychometricCurve();
    %NicePlot.exportFigToPDF(fullfile(videoOutDir,'PsychometricCurveComponent.pdf'), hFig, 300);

    visualizedWavelengths = [450 532 624 708];
    visualizedWavelengthsPSF = [450 480 532 550 624 708];
    
    theDisplay = generateTheDisplay();
    theMosaic = loadTheMosaic(rootPath, 16);
    [theScene, theRGBsettingsImage] = generateTheScene(rootPath, theDisplay, theMosaic.fov);
    
    hFig = generateDisplayFigComponent(theDisplay);
    NicePlot.exportFigToPDF(fullfile(videoOutDir,'/DisplayComponent.pdf'), hFig, 300);
    
    hFig = generateSceneFigComponent(theScene, theMosaic.fov, visualizedWavelengths);
    NicePlot.exportFigToPNG(fullfile(videoOutDir,'SceneComponent.png'), hFig, 300);
    
    hFig = generateRGBSettingsImageComponent(theRGBsettingsImage);
    NicePlot.exportFigToPNG(fullfile(videoOutDir,'RGBSettingsComponent.png'), hFig, 300);
 
    theOI = generateTheOpticalImage(theScene);
    theOI = oiCompute(theOI, theScene);
    
    hFig = generateOpticalImageFigComponent(theOI, theMosaic.fov, visualizedWavelengths);
    NicePlot.exportFigToPNG(fullfile(videoOutDir,'OpticalImageComponent.png'), hFig, 300);
    
    hFig = generatePSFFigComponent(theOI, visualizedWavelengthsPSF);
    NicePlot.exportFigToPNG(fullfile(videoOutDir, 'PSFsComponent.png'), hFig, 300);
    
    emPathLengthSecs = 0.45;
    theMosaic.integrationTime = 2.5/1000;
    fixationalEyeMovementsNum = round(emPathLengthSecs/theMosaic.integrationTime);
    nTrials = 1;
    theEMPath = generateTheEMPath(theMosaic, fixationalEyeMovementsNum, nTrials, videoOutDir);
 
    
    [theIsomerizations, thePhotocurrents, LMSfilters] = ...
        theMosaic.compute(theOI, 'currentFlag',true, 'emPath', theEMPath);
    
    hFig = generatePhotocurrentFilters(theMosaic, LMSfilters);
    NicePlot.exportFigToPNG(fullfile(videoOutDir,'PhotocurrentFilters.png'), hFig, 300);
    
    hFig = generateMosaicFigComponent(theMosaic);
    NicePlot.exportFigToPNG(fullfile(videoOutDir,'MosaicComponent.png'), hFig, 300);
    
    makeGif = true;
    hFig = generateResponseFigComponent(theMosaic, theIsomerizations, 'isomerizations', videoOutDir, makeGif);
    NicePlot.exportFigToPNG(fullfile(videoOutDir,'IsomerizationsComponent.png'), hFig, 300);
    
    hFig = generateResponseFigComponent(theMosaic, thePhotocurrents, 'photocurrents', videoOutDir, makeGif);
    NicePlot.exportFigToPNG(fullfile(videoOutDir,'PhotocurrentsComponent.png'), hFig, 300);
end

function hFig = generatePhotocurrentFilters(theMosaic, LMSfilters)
    hFig = figure(7); clf
    set(hFig, 'Position', [10 10 300 380], 'Color', [1 1 1]);
    ax = subplot('Position', [0.14 0.13 0.82 0.85]);
    timeAxis = (1:size(LMSfilters,1)) * theMosaic.integrationTime*1000;
    plot(ax,timeAxis, squeeze(LMSfilters(:,1)), 'r-', 'LineWidth', 2);
    hold on;
    plot(ax,timeAxis, squeeze(LMSfilters(:,2)), 'g-', 'Color', [0 0.75 0], 'LineWidth', 2);
    plot(ax,timeAxis, squeeze(LMSfilters(:,3)), 'b-', 'LineWidth', 2);
    hL = legend({'L', 'M', 'S'});
    set(ax, 'XLim', [0 400], 'YLim', [-0.02 0.16], 'YTick', [-0.02:0.02:0.2], 'YTickLabel', sprintf('%0.2f\n', -0.02:0.02:0.2));
    set(ax, 'FontSize', 14);
    xlabel(ax,'time (msec)', 'FontWeight', 'bold');
    grid(ax, 'on');
    box(ax, 'on');
end

function hFig = generateMosaicFigComponent(theMosaic)
    hFig = figure(4); clf
    set(hFig, 'Position', [10 10 426 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    
    theMosaic.visualizeGrid('axesHandle', ax, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'geometricArea', ...
        'backgroundColor', [1 1 1]);
    set(ax,'XColor', 'none', 'YColor', 'none');
    set(ax,'XTickLabels', {});
    set(ax,'YTickLabels', {});
    
end


function hFig = generatePSFFigComponent(theOI, visualizedWavelengths)
    micronsPerDegree = 300;
    hFig = visualizePSFfromOI(theOI, micronsPerDegree, ...
                'colormapToUse', gray(1024), ...
                'visualizedWavelengths', visualizedWavelengths, ...
                'rows', 2, 'cols', 3, ...
                'labelLastPSF', false, ...
                'displayWavelengthInTitle', ~false);
end

function hFig = generateResponseFigComponent(theMosaic, visualizedActivationPattern, signalName, videoOutDir, makeGif)

    % Get data for first trial
    nTrials = size(visualizedActivationPattern,1);
    visualizedActivationPatternTmp = squeeze(visualizedActivationPattern(1,:,:,:));
    emNum = size(visualizedActivationPatternTmp,3);
    
    % Compute signal range
    nonZeroConeIndices = find(theMosaic.pattern > 1);
    v = [];
    for emIndex = 1:emNum
        vv = squeeze(visualizedActivationPatternTmp(:,:,emIndex));
        v = cat(2,v,vv(nonZeroConeIndices));
    end
    signalRange = round(prctile(v(:), [5 95]));


    
    % Plot
    hFig = figure(5); clf

    if (makeGif)
        gifSizePixels = 96;
        minFigSize = 160;
        marginPixels = (minFigSize - gifSizePixels)/2;
        set(hFig, 'Position', [10 10 minFigSize/2 minFigSize/2], 'Color', [1 1 1], 'MenuBar', 'none', 'ToolBar', 'none', 'DockControls', 'off');
        ax = subplot('Position', [0.00 0.00 1.0 1.0]);
        gifFineName = fullfile(videoOutDir,sprintf('%sVideo.gif',signalName));
        showColorBar = false;
        labelColorBarTicks = false;
        backgroundColor = [0.2 0.2 0.2];
        blurGif = true;
        
    else
        set(hFig, 'Position', [10 10 400 420], 'Color', [1 1 1]);
        ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    
        mp4FileName = fullfile(videoOutDir,sprintf('%sVideo.mp4',signalName));
        videoOBJ = VideoWriter(mp4FileName, 'MPEG-4'); % H264 format
        videoOBJ.FrameRate = 30;
        videoOBJ.Quality = 100;
        videoOBJ.open();
        
        showColorBar = true;
        labelColorBarTicks = true;
        backgroundColor = [0 0 0];
    end
    
    activationLUT = gray(1024);
    titleForColorBar = ''; %'R*/5 msec';
    
    
    sampledHexMosaicXaxis = squeeze(theMosaic.patternSupport(1, :, 1));
    sampledHexMosaicYaxis = squeeze(theMosaic.patternSupport(:, 1, 2));
    xFullRange = max(sampledHexMosaicXaxis(:));
    yFullRange = max(sampledHexMosaicYaxis(:));
    
    
    
    frameIndex = 0;
    for iTrial = 1:nTrials
    for emIndex = 1:emNum
        frameIndex = frameIndex + 1;
        peakFrame = emNum/2;
        if (makeGif)
            ff = emIndex - peakFrame;
            sigma = peakFrame/3;
            fractionGauss = 2*exp(-0.5*(ff/sigma)^2);
            if (fractionGauss > 1)
                fractionGauss = 1;
            end
            minFraction = 0.55;
            maxFraction = 1.1;
            fraction = minFraction + fractionGauss;
            if (fraction > maxFraction)
                fraction = maxFraction;
            end

            xRange = xFullRange * fraction*[-1 1]/1e-6;
            yRange = yFullRange * fraction*[-1 1]/1e-6;
        else
            xRange = [];
            yRange = [];
        end
    
        theMosaic.renderActivationMap(ax, squeeze(visualizedActivationPattern(iTrial,:,:,emIndex)), ...
             'signalRange', signalRange, ...
             'visualizedConeAperture', 'geometricArea', ...
             'mapType', 'modulated disks', ...
             'showColorBar', showColorBar, ...
             'labelColorBarTicks', labelColorBarTicks, ...
             'xRange', xRange, ...
             'yRange', yRange, ...
             'showXLabel', false, ...
             'showYLabel', false, ...
             'titleForColorBar', titleForColorBar, ...
             'colorMap', activationLUT, ...
             'backgroundColor', backgroundColor);
        ylabel(ax, '');
        axis(ax, 'ij');
        xlim(ax, xRange*1e-6);
        ylim(ax, yRange*1e-6);
    
        set(ax,'XTickLabels', {});
        set(ax,'YTickLabels', {});
        
        drawnow;
        % Capture frame
        frame = getframe(hFig);

        if (makeGif)
            
            % Smooth with the smoothing kernel
            data = frame.cdata;
            rFactor = fractionGauss;
            gaussSigmaPixels = 0.9-0.2*rFactor^0.5; % 0.6*(2.0*rFactor^0.5);
            frame.cdata = data(marginPixels+(1:gifSizePixels), marginPixels+(1:gifSizePixels),:);
            data = frame.cdata;
            
            if (blurGif)
                 dataGray = squeeze(sum(data,3))/3;
   
                 dataGray = imgaussfilt(dataGray, gaussSigmaPixels);
                 maxNow = prctile(dataGray(:), 98-3*rFactor^1.0);
                 dataGray = 256 * dataGray / maxNow;
                 dataGray(dataGray>256) = 256;
 
                 for plane = 1:3
                     data(:,:,plane) = uint8(dataGray);
                 end
                 frame.cdata = data;
            end
             
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256, 'nodither');
            % Write to the GIF File 
              if frameIndex == 1 
                  imwrite(imind,cm,gifFineName,'gif', ...
                      'Loopcount',inf, ...
                      'DelayTime', 0/1000); 
              else 
                  imwrite(imind,cm,gifFineName,'gif', ...
                      'WriteMode','append', ...
                      'DelayTime', 0/1000); 
              end 
        else
            % Add video frame
            videoOBJ.writeVideo(frame);
        end
    end
    end
    
    if (makeGif)
        % view it
        web(gifFineName);
    else
        videoOBJ.close();
    end
    
end

function hFig = generateOpticalImageFigComponent(theOI, mosaicFOV, visualizedWavelengths)
    rgbImage = oiGet(theOI, 'rgb image');
    sampleSizeDegs = oiGet(theOI, 'wangular resolution');
    xSupport = 1:size(rgbImage,2);
    xSupport = xSupport * sampleSizeDegs;
    xSupport = xSupport - mean(xSupport);
    ySupport = 1:size(rgbImage,1);
    ySupport = ySupport * sampleSizeDegs;
    ySupport = ySupport - mean(ySupport);
    xyRange = [max(xSupport) max(ySupport)]*0.8*2.0;
    
    wave = oiGet(theOI, 'wave');
    photons = oiGet(theOI, 'photons');
    maxPhotons = max(photons(:));
    minPhotons = min(photons(:));
    
    hFig = figure(3); clf;
    set(hFig, 'Position', [10 10 500 500], 'Color', [1 1 1]);

    renderedImagesNum = 1+numel(visualizedWavelengths);
    dx = 0.64/renderedImagesNum;
    dy = 0.68/renderedImagesNum;
    for k = 1:renderedImagesNum
        if (k < renderedImagesNum)
            targetWavelength = visualizedWavelengths(k);
            [~,idx] = min(abs(wave-targetWavelength));
            thePhotons = squeeze(photons(:,:,idx));
            imageData = (thePhotons-minPhotons)/(maxPhotons-minPhotons);
        else
            targetWavelength = [];
            imageData = rgbImage;
        end
        
        ax = axes('Position', [0.01+(k-1)*dx 0.51-(k-1)*dy 0.47 0.53]);
        renderImage(ax,xSupport, ySupport, imageData, mosaicFOV, xyRange, targetWavelength);

        drawnow
    end
end

function hFig = generateRGBSettingsImageComponent(RGBsettingsImage)
    hFig = figure(3); clf;
    set(hFig, 'Position', [10 10 400 500], 'Color', [1 1 1]);
    xx = 1:size(RGBsettingsImage,2);
    yy = 1:size(RGBsettingsImage,1);
    renderedImagesNum = 3;
    dx = 0.12/renderedImagesNum;
    dy = 0.68/renderedImagesNum;
    for k = 1:renderedImagesNum
        imageData = RGBsettingsImage;
        switch (k)
            case 1
                imageData(:,:,2) = 0;
                imageData(:,:,3) = 0;
            case 2
                imageData(:,:,1) = 0;
                imageData(:,:,3) = 0;
            case 3
                imageData(:,:,1) = 0;
                imageData(:,:,2) = 0;
        end
        
        ax = axes('Position', [0.01+(k-1)*dx 0.40-(k-1)*dy 0.9 0.65]);
        imagesc(ax, xx, yy, imageData, [0 1]);
        set(gca, 'XTick', [], 'YTick', []);
        axis image
    end
end

function hFig = generateSceneFigComponent(theScene, mosaicFOV, visualizedWavelengths)

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
    
    wave = sceneGet(theScene, 'wave');
    photons = sceneGet(theScene, 'photons');
    maxPhotons = max(photons(:));
    minPhotons = min(photons(:));
    
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 500 500], 'Color', [1 1 1]);
    renderedImagesNum = 1+numel(visualizedWavelengths);
    dx = 0.64/renderedImagesNum;
    dy = 0.68/renderedImagesNum;
    for k = 1:renderedImagesNum
        if (k < renderedImagesNum)
            targetWavelength = visualizedWavelengths(k);
            [~,idx] = min(abs(wave-targetWavelength));
            thePhotons = squeeze(photons(:,:,idx));
            imageData = (thePhotons-minPhotons)/(maxPhotons-minPhotons);
        else
            targetWavelength = [];
            imageData = sceneGet(theScene, 'rgb image');
        end
        
        ax = axes('Position', [0.01+(k-1)*dx 0.51-(k-1)*dy 0.47 0.53]);
        renderImage(ax,xSupport, ySupport, imageData, mosaicFOV, sceneFOV, targetWavelength);

        drawnow
    end
end

function renderImage(ax,xSupport, ySupport, imageData, mosaicFOV, xyRange, targetWavelength)
    imagesc(ax, xSupport, ySupport, imageData, [0 1]);
    hold(ax, 'on');
    x = 0.5*mosaicFOV(1) * [-1 1 1 -1 -1];
    y = 0.5*mosaicFOV(2) * [-1 -1 1 1 -1];
    plot(ax,x,y,'k--', 'LineWidth', 1.0);
    
    axis(ax, 'image');% axis(ax, 'xy');
    tickInterval = 0.2;
    xLim = 0.75*0.4*xyRange(1)*[-1 1];
    yLim = 0.75*0.5*xyRange(2)*[-1 1];
    outline.x = [xLim(1) xLim(1) xLim(2) xLim(2) xLim(1)];
    outline.y = [yLim(1) yLim(2) yLim(2) yLim(1) yLim(1)];
    plot(outline.x, outline.y, 'k-', 'LineWidth', 1.0);
    hold(ax, 'off');
    set(ax, 'XLim', xLim, 'YLim', yLim, ...
         'XTick', -15:tickInterval:15, 'YTick', -15:tickInterval:15, ...
         'XTickLabels', {}, 'YTickLabels', {});
    grid(ax, 'off'); box(ax, 'on');
    set(ax, 'FontSize', 12);
    
    if (~isempty(targetWavelength))
        xT = xLim(1) + 0.005;
        yT = yLim(1) + 0.025;
        t = text(xT,yT,sprintf('\\lambda = %2.0f nm', targetWavelength)); 
        set(t, 'Color', [1 1 1], 'FontSize', 20);
        colormap(ax, gray(1024));
    end
    
end

function hFig = generateDisplayFigComponent(display)
    hFig = figure(1); clf
    set(hFig, 'Position', [10 10 426 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.15 0.96 0.85]);
    
    wave = [display.wave(1)];
    red = [0]; green = [0]; blue = [0];
    for k = 1:numel(display.wave)
        wave = [wave display.wave(k)];
        red = [red display.spd(k,1)];
        green = [green display.spd(k,2)];
        blue = [blue display.spd(k,3)];
        if (k < numel(display.wave))
            wave = [wave display.wave(k+1)];
            red = [red display.spd(k,1)];
            green = [green display.spd(k,2)];
            blue = [blue display.spd(k,3)];
        end
    end
    
    area(ax, wave, red, 'EdgeColor',[1 0.0 0.0], 'LineWidth', 1.5, 'FaceColor',[1 0.5 0.5],'FaceAlpha',.5,'EdgeAlpha',.8);
    hold(ax, 'on');
    area(ax, wave, green, 'EdgeColor',[0 1.0 0.0], 'LineWidth', 1.5, 'FaceColor',[0.5 1.0 0.5],'FaceAlpha',.5,'EdgeAlpha',.8);
    area(ax, wave, blue, 'EdgeColor',[0 0.0 1.0], 'LineWidth', 1.5, 'FaceColor',[0.5 0.5 1],'FaceAlpha',.5,'EdgeAlpha',.8);
    grid off;
    set(gca, 'XLim', [display.wave(1) display.wave(end)]); 
    set(gca, 'Color', [1 1 1], 'XTick', [400:50:800], 'YTick', [], 'YColor', 'none', 'XColor', [0 0 0]);
    xlabel('wavelength, \lambda (nm)', 'Color', [0 0 0]);
    set(gca, 'FontSize', 20);
    box off
    
end

function theEMPath = generateTheEMPath(theMosaic, eyeMovementsNum, nTrials, videoOutDir)
    fixEMobj = fixationalEM();
    randomSeed = 7;  % 857;
    fixEMobj.computeForConeMosaic(theMosaic, eyeMovementsNum, ...
        'nTrials', nTrials, ...
        'rSeed', randomSeed );
    theEMPath = fixEMobj.emPos;

    centerPath = true;
    if (centerPath)
        xMedian = median(squeeze(theEMPath(1,:,1)));
        yMedian = median(squeeze(theEMPath(1,:,2)));
        theEMPath(:,:,1) = theEMPath(:,:,1) - xMedian;
        theEMPath(:,:,2) = theEMPath(:,:,2) - yMedian;
        % We cannot set the emPos, so the visualization below shows the 
        % uncentered path
        %fixEMobj.emPos = theEMPath;
    end
    
    if (1==1)
    fixEMobj.generateEMandMosaicComboVideo(...
            fixEMobj, theMosaic, ...
            'visualizedFOVdegs', 0.8*theMosaic.fov(1), ...
            'showMovingMosaicOnSeparateSubFig', true, ...
            'displaycrosshairs', ~true, ...
            'videoFileName', fullfile(videoOutDir, 'fixationalEMVideoNoCrossHairs.mp4'));
    end
        
end

function theOI = generateTheOpticalImage(theScene)
    opticsModel = 'ThibosDefaultSubject3MMPupil';  % 'Geisler'; % '
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


function theMosaic = loadTheMosaic(rootPath, mosaicFOVDegs)
    mosaicFile = sprintf('coneMosaic%dCPD.mat', mosaicFOVDegs);
    filename = fullfile(rootPath,mosaicFile);
    load(filename, 'theData');
    theMosaic = theData;
end

function [theScene, sensorImageRGB] = generateTheScene(rootPath,theDisplay, theMosaicFOV)
%   
    rgbSettingsFile = 'PilotC4M4-RGB.mat';  % matte
    rgbSettingsFile = 'PilotC1M4-RGB.mat';  % shiny
    filename = fullfile(rootPath,rgbSettingsFile);
    load(filename, 'sensorImageRGB');
    
    %% Convert to ISETbio scene
    meanLuminance = [];
    theScene = sceneFromFile(sensorImageRGB, 'rgb', meanLuminance, theDisplay);
        
    meanLuminance = sceneGet(theScene, 'mean luminance')
    theScene = sceneAdjustLuminance(theScene, 100)
    meanLuminance = sceneGet(theScene, 'mean luminance')
    pause
    % Scale the scene to be a little larger than the mosaic
    theScene = sceneSet(theScene, 'h fov', theMosaicFOV(1)*0.9);
end

function theDisplay = generateTheDisplay()
    backgroundParams.backgroundxyY = [0.27 0.30 49.8]';
    backgroundParams.monitorFile = 'CRT-MODEL';
    backgroundParams.leakageLum = 1.0;
    backgroundParams.lumFactor = 1;
    
    theDisplay = displayCreate(backgroundParams.monitorFile);
end