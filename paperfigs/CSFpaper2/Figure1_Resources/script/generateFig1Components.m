function generateFig1Components
    % This is the function used to make the components of Fig1  for the paper2 second submission.
    
    dataFileName = 'dynamicData.mat';
    
    redoComputations = ~true;
    if (redoComputations)
        
        [rootPath,~] = fileparts(which(mfilename));
        rootPath = strrep(rootPath, 'script', 'isetbio_resources');

        generateDisplayFig(rootPath);

        mosaicFOV = 0.6;
        theConeMosaic = generateMosaicFig(rootPath, mosaicFOV);

        % The scene
        sceneFOV = 1.0;
        visualizedSceneFraction = [0.6*1.1 0.6];
        theScene = generateSceneFig(rootPath, sceneFOV, visualizedSceneFraction, mosaicFOV);


        opticsModel = 'ThibosDefaultSubject3MMPupil'; pupilDiamMm = 3.0;
        theOI = generateOpticsImageFig(opticsModel, pupilDiamMm, sceneFOV);

        theOI = oiCompute(theOI, theScene);
        generateOpticalImageFig(theOI, sceneFOV*visualizedSceneFraction, mosaicFOV);

        % Static isomerizations
        theIsomerizations = theConeMosaic.compute(theOI, 'currentFlag',false);
        generateIsomerizationsImageFig(theConeMosaic, theIsomerizations);
        
        
        % Generate fixational eye movements for 1 trial and for the entire simulation time
        eyeMovementsNum = 120;
        [theEMPaths, theFixationalEMobject] = theConeMosaic.emGenSequence(eyeMovementsNum, 'nTrials', 10);
        
        % Select a nice emPath (9 is a good one)
        selectedEMPathIndex = 9; % selectEMPath(theEMPaths, theConeMosaic);
        
        % Generate a number of identical emPaths(selectedEMPathIndex)
        theEMPaths = squeeze(theEMPaths(selectedEMPathIndex,:,:));
        nTrialsNum = 4;
        theEMPaths = repmat(reshape(theEMPaths, [1 eyeMovementsNum 2]), [nTrialsNum 1 1]);

        % Dynamic isomerizations
        % Compute mosaic responses to the same emPath with noise ON
        theConeMosaic.noiseFlag = 'random';
        theConeMosaic.os.noiseFlag = 'random';
        
        % Compute mosaic responses to the same repeated emPath with noise ON
        [coneExcitations, photoCurrents] = theConeMosaic.compute(theOI, ...
            'emPath', theEMPaths, 'currentFlag', true);
        
        % Compute mosaic responses to the same emPath with noise OFF
        theConeMosaic.noiseFlag = 'none';
        theConeMosaic.os.noiseFlag = 'none';
        [coneExcitationsMean, photoCurrentsMean] = theConeMosaic.compute(theOI, ...
            'emPath', theEMPaths(1:2,:,:), 'currentFlag', true);
        
        % Save the data
        save(dataFileName, 'theConeMosaic', 'theFixationalEMobject', 'theEMPaths',  'coneExcitations', 'photoCurrents', 'coneExcitationsMean', 'photoCurrentsMean', '-v7.3');
    else
        % Load the data
        load(dataFileName, 'theConeMosaic', 'theFixationalEMobject', 'theEMPaths',  'coneExcitations', 'photoCurrents', 'coneExcitationsMean', 'photoCurrentsMean');
        
        % Visualize the first EMPath on the cone mosaic using a
        % color-varying line
        visualizedTrialIndex = 1;
        visualizedConePosMicrons = [-84.74 -43.77]; % [-75.1 -55.16]; % [37.31 21.32];
        visualizedTimeRange = [0 590];
        visualizedConePosDegs = visualizedConePosMicrons / theConeMosaic.micronsPerDegree;
        [coneIndices, conePositions, coneTypes] = theConeMosaic.indicesForConesAtPositions(visualizedConePosDegs);
        visualizedConeIndex = coneIndices(1);
        visualizedConePosDegs = conePositions(1,:);
        visualizedConePosMicrons = visualizedConePosDegs * theConeMosaic.micronsPerDegree;
        theEMPathsMeters = squeeze(theEMPaths(visualizedTrialIndex,:,:)) * theConeMosaic.patternSampleSize(1);
        
        tBins = size(theEMPathsMeters,1);
        cMap = brewermap(tBins, 'Spectral');
        generateMosaicWithEMPathFig(theConeMosaic, theEMPathsMeters, tBins, cMap);
        
        % Visualize the first cone response mosaic instance at specific times together with the EMPath
        depictedTimes = 0:50:550;
        timeAxis = (0:1:(tBins-1))*theConeMosaic.integrationTime*1000;
        visualizedTrialConeExcitationResponseXYT = squeeze(coneExcitations(visualizedTrialIndex,:,:,:));
        
        videoFileName = 'response.mp4';
        videoOBJ = VideoWriter(videoFileName, 'MPEG-4'); % H264 format
        videoOBJ.FrameRate = 30;
        videoOBJ.Quality = 100;
        videoOBJ.open();
        
        for k = 1:numel(depictedTimes)
            [~,tBin] = min(abs((timeAxis-depictedTimes(k))));
            hFig = generateMosaicResponseWithEMPathFig(100+k,theConeMosaic, theEMPathsMeters, visualizedTrialConeExcitationResponseXYT, visualizedConePosMicrons, tBin, cMap, timeAxis(tBin));
            NicePlot.exportFigToPDF(sprintf('componentFigs/MosaicResponseWithEMPathComponent_%dmsec.pdf',timeAxis(tBin)), hFig, 300);
            videoOBJ.writeVideo(getframe(hFig));
        end
        videoOBJ.close();
        fprintf('File saved in %s\n', videoFileName);
        
        
        % % Collect all responses for the selected cone
        trialsNum = size(coneExcitations,1);
        visualizedConeExcitationResponses = zeros(trialsNum, tBins);
        visualizedPhotocurrentResponses = zeros(trialsNum, tBins);
        visualizedMeanConeExcitationResponse = zeros(1, tBins);
        visualizedMeanPhotocurrentResponse = zeros(1, tBins);
        
        for trialIndex = 1:trialsNum
            tmp = squeeze(coneExcitations(trialIndex,:,:,:));
            tmp = reshape(tmp, [size(tmp,1)*size(tmp,2) size(tmp,3)]);
            visualizedConeExcitationResponses(trialIndex,:) = tmp(visualizedConeIndex,:);
            
            tmp = squeeze(photoCurrents(trialIndex,:,:,:));
            tmp = reshape(tmp, [size(tmp,1)*size(tmp,2) size(tmp,3)]);
            visualizedPhotocurrentResponses(trialIndex,:) = tmp(visualizedConeIndex,:);
            
            tmp = squeeze(coneExcitationsMean(1,:,:,:));
            tmp = reshape(tmp, [size(tmp,1)*size(tmp,2) size(tmp,3)]);
            visualizedMeanConeExcitationResponse = tmp(visualizedConeIndex,:);
            
            tmp = squeeze(photoCurrentsMean(1,:,:,:));
            tmp = reshape(tmp, [size(tmp,1)*size(tmp,2) size(tmp,3)]);
            visualizedMeanPhotocurrentResponse = tmp(visualizedConeIndex,:);
        end
        
        % Visualize the responses of the selected cone
        figNo = 1000;
        generateSingleConeResponsePlots(figNo, timeAxis, visualizedTimeRange, ...
            visualizedMeanConeExcitationResponse, visualizedConeExcitationResponses, ...
            visualizedMeanPhotocurrentResponse, visualizedPhotocurrentResponses, cMap);
        
    end       
end      
     
function generateSingleConeResponsePlots(figNo, timeAxis, visualizedTimeRange, ...
    visualizedMeanConeExcitationResponse, visualizedConeExcitationResponses, ...
    visualizedMeanPhotocurrentResponse, visualizedPhotocurrentResponses, cMap)

    hFig = figure(figNo); clf
    set(hFig, 'Position', [10 10 426 700], 'Color', [1 1 1]);
    
    % Render cone excitation responses at the top
    XLim = visualizedTimeRange;
    XTick = 0:100:600;
    YRange = max(visualizedMeanConeExcitationResponse)-min(visualizedMeanConeExcitationResponse);
    YLim = min(visualizedMeanConeExcitationResponse) + YRange*[-0.25 1.25];
    YTick = 0:10:100;
    ax = subplot('Position', [0.14 0.55 0.84 0.44]);
    plotMeanAndSingleResponseInstances(ax, timeAxis, visualizedMeanConeExcitationResponse, visualizedConeExcitationResponses, ...
        XLim, YLim, XTick, YTick, 'cone excitations (R*/c/5 msec)', cMap, ~true);


    % Render pCurrent responses at the bottom
    YRange = max(visualizedMeanPhotocurrentResponse)-min(visualizedMeanPhotocurrentResponse);
    YLim = min(visualizedMeanPhotocurrentResponse) + YRange*[-0.25 1.25];
    YTick = -80:10:0;
    ax = subplot('Position', [0.14 0.07 0.84 0.44]);
    plotMeanAndSingleResponseInstances(ax, timeAxis, visualizedMeanPhotocurrentResponse, visualizedPhotocurrentResponses, ...
        XLim, YLim, XTick, YTick, 'photocurrent (pAmps)', cMap, true);
        
end


function plotMeanAndSingleResponseInstances(ax, timeAxis, meanResponse, responseInstances, ...
            XLim, YLim, XTick, YTick, yAxisLabel, cMap, displayXLabel)
        
    tBins = numel(timeAxis);
    plot(timeAxis, responseInstances, 'k-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); hold on
    plot(ax, timeAxis, meanResponse, 'k-', 'LineWidth', 5);
    % Overlay the trace in color
    for k = 1:tBins-1
       plot([timeAxis(k) timeAxis(k+1)], [meanResponse(k) meanResponse(k+1)], ...
          '-', 'Color', squeeze(cMap(k,:)), 'LineWidth', 2);
    end
    set(gca, 'XLim', XLim, 'XTick', XTick, 'YLim', YLim, 'YTick', YTick, 'FontSize', 16);
    grid on; box off;
    if (displayXLabel)
        xlabel('\it time (msec)');
    else
        set(gca, 'XTickLabel', {});
    end
    ylabel(sprintf('\\it %s', yAxisLabel));
end


function hFig = generateMosaicResponseWithEMPathFig(figNo, theConeMosaic, theEMPathsMeters, visualizedTrialConeExcitationResponseXYT, visualizedConePosMicrons, tBin, cMap, timeMsec)
    hFig = figure(figNo); clf
    set(hFig, 'Position', [10 10 426 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    
    theConeMosaic.renderActivationMap(ax, squeeze(visualizedTrialConeExcitationResponseXYT(:,:,tBin)),...
                'outlineConesAlongHorizontalMeridian', ~true, ...
                'visualizedConeAperture', 'geometricArea');
    hold on;
    plot(ax,visualizedConePosMicrons(1)*1e-6, visualizedConePosMicrons(2)*1e-6, 'cs', 'LineWidth', 2.0, 'MarkerSize', 16);
    plotCrossHairs(ax, theEMPathsMeters, tBin, cMap);
    textPosMicrons = [40 -80];
    text(ax, textPosMicrons(1)*1e-6, textPosMicrons(2)*1e-6, sprintf('%2.0f msec', timeMsec), 'Color', 'y', 'FontSize', 20);
    hold off;
    set(ax,'XColor', 'none', 'YColor', 'none');
    set(ax,'XTickLabels', {});
    set(ax,'YTickLabels', {});
 end



function selectedEMPathIndex = selectEMPath(theEMPaths, theConeMosaic)
        
    eyeMovementsNum = size(theEMPaths,2);
    timeAxis = (0:1:(eyeMovementsNum-1))*theConeMosaic.integrationTime*1000;

    depictedTimes = [0 50 100 150 200 250 300 350 400];
    cMap = brewermap(eyeMovementsNum, 'Spectral');
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 5, ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.07, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.05);

    for emPathIndex = 1:size(theEMPaths,1)
        theEMPath = squeeze(theEMPaths(emPathIndex,:,:));
        theEMPathMeters = theEMPath * theConeMosaic.patternSampleSize(1);
        
        % Show the 2D activation for a single response instance at different time epochs
        hFig = figure(emPathIndex); clf;
        set(hFig, 'Position', [10 + emPathIndex*10, 10 + emPathIndex*30, 2400 850]);
        for k = 1:numel(depictedTimes)
            row = floor((k-1)/5)+1;
            col = mod(k-1,5)+1;
            ax = subplot('Position', subplotPosVectors(row,col).v);
            [~,tBin] = min(abs((timeAxis-depictedTimes(k))));
            theConeMosaic.visualizeGrid('axesHandle', ax, ...
                'visualizedConeAperture', 'geometricArea');
            hold on;
            plotCrossHairs(ax, theEMPathMeters, tBin, cMap);
            hold off;
            title(sprintf('%d ms', timeAxis(tBin)));
            drawnow;
        end
    end
    
    selectedEMPathIndex = input('Enter desired emPath: ');
end

function plotCrossHairs(ax, theEMPath, tBin, cMap)
    % Plot the cross-hairs at (0,0)
    plot(ax, 100*1e-6*[-1 1], [0 0], 'w-', 'LineWidth', 1.5);
    plot(ax, [0 0], 100*1e-6*[-1 1], 'w-', 'LineWidth', 1.5);
    
    % Plot the trace in black
    plot(ax, theEMPath(1:tBin, 1), -theEMPath(1:tBin, 2), 'k-', 'LineWidth', 5);
    
    % Overlay the trace in color
    for k = 1:tBin-1
        plot(ax, [theEMPath(k, 1) theEMPath(k+1, 1)], ...
            -[theEMPath(k, 2) theEMPath(k+1, 2)] , '-', 'Color', squeeze(cMap(k,:)), 'LineWidth', 2);
    end
end


function generateMosaicWithEMPathFig(theConeMosaic, theEMPathsMeters, tBins, cMap)
    hFig = figure(1004); clf
    set(hFig, 'Position', [10 10 426 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    
    theConeMosaic.visualizeGrid('axesHandle', ax, ...
        'visualizedConeAperture', 'geometricArea', ...
        'noXaxisLabel', true, 'noYaxisLabel', true);
    hold on;
    plotCrossHairs(ax, theEMPathsMeters, tBins, cMap);
    hold off;
    set(ax,'XColor', 'none', 'YColor', 'none');
    set(ax,'XTickLabels', {});
    set(ax,'YTickLabels', {});
    NicePlot.exportFigToPDF('componentFigs/MosaicWithEMPathComponent.pdf', hFig, 300);
end




function generateDisplayFig(rootPath)
    load(fullfile(rootPath, 'display.mat'), 'display');
    hFig = figure(4); clf
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
    grid on;
    set(gca, 'XLim', [display.wave(1) display.wave(end)]); 
    set(gca, 'XTick', [400:50:800], 'YTick', [], 'YColor', 'none');
    xlabel('wavelength, \lambda (nm)');
    set(gca, 'FontSize', 20);
    box off
    pwd
    NicePlot.exportFigToPDF('componentFigs/DisplayComponent.pdf', hFig, 300);
end

function generateIsomerizationsImageFig(theConeMosaic, theIsomerizations)

    visualizedActivationPattern = theIsomerizations;
    signalRange = [0 max(visualizedActivationPattern(:))];
    
    hFig = figure(3); clf
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
     set(ax,'XTickLabels', {});
     set(ax,'YTickLabels', {});
     NicePlot.exportFigToPDF('componentFigs/MeanIsomerizationsComponent.pdf', hFig, 300);
end


function generateOpticalImageFig(theOI, xyRange, mosaicFOV)
    rgbImage = oiGet(theOI, 'rgb image');
    sampleSizeDegs = oiGet(theOI, 'wangular resolution');
    xSupport = 1:size(rgbImage,2);
    xSupport = xSupport * sampleSizeDegs;
    xSupport = xSupport - mean(xSupport);
    ySupport = 1:size(rgbImage,2);
    ySupport = ySupport * sampleSizeDegs;
    ySupport = ySupport - mean(ySupport);
    
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 468 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    image(ax,xSupport , ySupport, rgbImage);
    hold(ax, 'on');
    x = 0.5*mosaicFOV * [-1 1 1 -1 -1];
    y = 0.5*mosaicFOV * [-1 -1 1 1 -1];
    plot(ax, x,y,'k-', 'LineWidth', 1.0);
    hold(ax, 'off');
    axis(ax, 'xy'); axis(ax, 'image');
    
    tickInterval = 0.2;
    set(ax, 'XLim', xyRange(1)/2*[-1 1], 'YLim', xyRange(2)/2*[-1 1], ...
         'XTick', -15:tickInterval:15, 'YTick', -15:tickInterval:15, ...
         'YTickLabels', {}, 'YTickLabels', {});
    grid on; box on;
    set(ax, 'FontSize', 12);
    %xlabel('space (degs)');
    %ylabel('space (degs)');
    NicePlot.exportFigToPNG('componentFigs/OpticalImageComponent.png', hFig, 300);
    
    
    
    wave = oiGet(theOI, 'wave');
    photons = oiGet(theOI, 'photons');
    maxPhotons = max(photons(:));
    minPhotons = min(photons(:));
    visualizedWaves = [450 480 510 530 550 580 610 640 670];
    for k = 1:numel(visualizedWaves)
        [~,idx] = min(abs(wave-visualizedWaves(k)));
        thePhotons = squeeze(photons(:,:,idx));
        hFig = figure(2999+k); clf;
        set(hFig, 'Position', [10 10 468 420], 'Color', [1 1 1]);
        ax = subplot('Position', [0.02 0.02 0.96 0.96]);
        imagesc(ax, xSupport, ySupport, (thePhotons-minPhotons)/(maxPhotons-minPhotons), [0 1]);
        hold(ax, 'on');
        x = 0.5*mosaicFOV * [-1 1 1 -1 -1];
        y = 0.5*mosaicFOV * [-1 -1 1 1 -1];
        plot(ax,x,y,'k-', 'LineWidth', 1.0);
        hold(ax, 'off');
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'XLim', xyRange(1)/2*[-1 1], 'YLim', xyRange(2)/2*[-1 1], ...
         'XTick', -15:tickInterval:15, 'YTick', -15:tickInterval:15, ...
         'XTickLabels', {}, 'YTickLabels', {});
        colormap(gray(1024));
        NicePlot.exportFigToPNG(sprintf('componentFigs/PhotonsOIComponent%d.png', wave(idx)), hFig, 300);
    end
    
    
    micronsPerDegree = 300;
    visualizePSFfromOI(theOI, micronsPerDegree, ...
                'colormapToUse', gray(1024), ...
                'visualizedWavelengths', [450 480 510 530 550 580 610 640 670], ...
                'rows', 3, 'cols', 3, ...
                'labelLastPSF', false, ...
                'displayWavelengthInTitle', ~false);
            
   
end



function theScene = generateSceneFig(rootPath, sceneFOV, visualizedSceneFraction, mosaicFOV)
    load(fullfile(rootPath, 'scenes.mat'));
    theScene = scenesList{1};
    
    theScene = sceneSet(theScene, 'h fov', sceneFOV);
    sceneGet(theScene, 'wAngular')
    
    meanLuminance = 34;
    theScene = sceneAdjustLuminance(theScene, meanLuminance);
    meanLuminance = sceneGet(theScene, 'mean luminance');
    
    save(fullfile('componentFigs','scene.mat'), 'theScene');
    
    hFig = figure(1); clf
    set(hFig, 'Position', [10 10 468 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    sceneSize = sceneGet(theScene, 'size');
    rows = sceneSize(1); cols = sceneSize(2);
    sceneGet(theScene, 'width', 'm')
    xSupport = (1:cols);
    xSupport = xSupport - mean(xSupport);
    xSupport = 0.5*xSupport / max(abs(xSupport));
    xSupport = xSupport * sceneGet(theScene, 'wAngular');
    ySupport = (1:rows);
    ySupport = ySupport - mean(ySupport);
    ySupport = 0.5*ySupport / max(abs(ySupport));
    ySupport = ySupport * sceneGet(theScene, 'hAngular');
    
    imagesc(ax, xSupport, ySupport, sceneGet(theScene, 'rgb image'));
    hold(ax, 'on');
    x = 0.5*mosaicFOV * [-1 1 1 -1 -1];
    y = 0.5*mosaicFOV * [-1 -1 1 1 -1];
    plot(ax,x,y,'k-', 'LineWidth', 1.0);
    hold(ax, 'off');
    
    axis(ax, 'image'); axis(ax, 'xy');
    tickInterval = 0.2;
    xyRange = sceneFOV * visualizedSceneFraction;
    set(ax, 'XLim', xyRange(1)/2*[-1 1], 'YLim', xyRange(2)/2*[-1 1], ...
         'XTick', -15:tickInterval:15, 'YTick', -15:tickInterval:15, ...
         'XTickLabels', {}, 'YTickLabels', {});
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontSize', 12);
    %xlabel('space (degs)');
    %ylabel('space (degs)');
    NicePlot.exportFigToPNG('componentFigs/SceneComponent.png', hFig, 300);
    
    
    wave = sceneGet(theScene, 'wave');
    photons = sceneGet(theScene, 'photons');
    maxPhotons = max(photons(:));
    minPhotons = min(photons(:));
    visualizedWaves = [450 480 510 530 550 580 610 640 670];
    for k = 1:numel(visualizedWaves)
        [~,idx] = min(abs(wave-visualizedWaves(k)));
        thePhotons = squeeze(photons(:,:,idx));
        hFig = figure(1999+k); clf;
        set(hFig, 'Position', [10 10 468 420], 'Color', [1 1 1]);
        ax = subplot('Position', [0.02 0.02 0.96 0.96]);
        imagesc(ax, xSupport, ySupport, (thePhotons-minPhotons)/(maxPhotons-minPhotons), [0 1]);
        hold(ax, 'on');
        x = 0.5*mosaicFOV * [-1 1 1 -1 -1];
        y = 0.5*mosaicFOV * [-1 -1 1 1 -1];
        plot(ax,x,y,'k-', 'LineWidth', 1.0);
        hold(ax, 'off');
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'XLim', xyRange(1)/2*[-1 1], 'YLim', xyRange(2)/2*[-1 1], ...
         'XTick', -15:tickInterval:15, 'YTick', -15:tickInterval:15, ...
         'XTickLabels', {}, 'YTickLabels', {});
        colormap(gray(1024));
        NicePlot.exportFigToPNG(sprintf('componentFigs/PhotonsSceneComponent%d.png', wave(idx)), hFig, 300);
    end
    
    
end

function theConeMosaic = generateMosaicFig(rootPath, mosaicFOV)
    load(fullfile(rootPath, sprintf('coneMosaic_%1.2fdegFOV.mat', mosaicFOV)), 'theConeMosaic');
    
    hFig = figure(4); clf
    set(hFig, 'Position', [10 10 426 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    
    theConeMosaic.visualizeGrid('axesHandle', ax, ...
        'visualizedConeAperture', 'geometricArea');
    set(ax,'XColor', 'none', 'YColor', 'none');
    set(ax,'XTickLabels', {});
    set(ax,'YTickLabels', {});
    NicePlot.exportFigToPDF('componentFigs/MosaicComponent.pdf', hFig, 300);
end


function theOI = generateOpticsImageFig(opticsModel, pupilDiamMm, horizontalFOV)
    oiParams = struct(...
                'opticsModel', opticsModel, ...
                'wavefrontSpatialSamples', 301, ...
                'pupilDiamMm', pupilDiamMm, ...
                'umPerDegree', 300);
    theOI = oiWithCustomOptics(oiParams.opticsModel, oiParams.wavefrontSpatialSamples, oiParams.pupilDiamMm, oiParams.umPerDegree);   
    
    % Set the FOV
    theOI = oiSet(theOI,'h fov',horizontalFOV);

    % Set the fNumber
    focalLength = oiGet(theOI,'distance');
    desiredFNumber = focalLength/(oiParams.pupilDiamMm/1000);
    theOI  = oiSet(theOI ,'optics fnumber',desiredFNumber);  
end


