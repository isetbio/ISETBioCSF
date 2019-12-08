function generateFig1Components
    % This is the function used to make the components of Fig1  for the paper2 second submission.
    nTrialsNum = 2;
    dataFileName = sprintf('dynamicDataTrials%1.0f.mat', nTrialsNum);
    
    redoComputations = ~true;
    if (redoComputations)
        
        [rootPath,~] = fileparts(which(mfilename));
        rootPath = strrep(rootPath, 'script', 'isetbio_resources');

        generateDisplayFig(rootPath);

        mosaicFOV = 0.6;
        theConeMosaic = generateMosaicFig(rootPath, mosaicFOV, []);

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
        theEMPaths = repmat(reshape(theEMPaths, [1 eyeMovementsNum 2]), [nTrialsNum 1 1]);

        % Dynamic isomerizations
        % Compute mosaic responses to the same emPath with noise ON
        theConeMosaic.noiseFlag = 'random';
        theConeMosaic.os.noiseFlag = 'random';
        
        % Compute mosaic responses to the same repeated emPath with noise ON
        [coneExcitations, photoCurrents, pCurrentImpulseResponses] = theConeMosaic.compute(theOI, ...
            'emPath', theEMPaths, 'currentFlag', true);
        
        % Compute mosaic responses to the same emPath with noise OFF
        theConeMosaic.noiseFlag = 'none';
        theConeMosaic.os.noiseFlag = 'none';
        [coneExcitationsMean, photoCurrentsMean] = theConeMosaic.compute(theOI, ...
            'emPath', theEMPaths(1:2,:,:), 'currentFlag', true);
        
        % Save the data
        save(dataFileName, 'theConeMosaic', 'theFixationalEMobject', 'theEMPaths',  'pCurrentImpulseResponses', 'coneExcitations', 'photoCurrents', 'coneExcitationsMean', 'photoCurrentsMean', '-v7.3');
    else
        % Load the data
        load(dataFileName, 'theConeMosaic', 'theFixationalEMobject', 'theEMPaths',  'pCurrentImpulseResponses', 'coneExcitations', 'photoCurrents', 'coneExcitationsMean', 'photoCurrentsMean');
        
        generateMosaicFig([], [], theConeMosaic);
        
        % Visualize the first EMPath on the cone mosaic using a
        % color-varying line
        visualizedTrialIndex = 1;
        visualizedTimeRange = [0 590];
        
        
        % Select 3 cones to visualize their responses
        visualizedConePosMicrons = [...
            -85 -46; ...
            -84.74 -43.77; ...
            -88 -43];
        
        
        % Find the corresponding cone indices and types
        for coneIdx = 1:3
            visualizedConePosDegs = squeeze(visualizedConePosMicrons(coneIdx,:)) / theConeMosaic.micronsPerDegree;
            [coneIndices, conePositions, coneTypes] = theConeMosaic.indicesForConesAtPositions(visualizedConePosDegs);
            visualizedConeIndex(coneIdx) = coneIndices(1);
            visualizedConeType(coneIdx) = coneTypes(1);
            visualizedConePosDegs = conePositions(1,:);
            visualizedConePosMicrons(coneIdx,:) = visualizedConePosDegs * theConeMosaic.micronsPerDegree;
        end

        
        theEMPathsMeters = squeeze(theEMPaths(visualizedTrialIndex,:,:)) * theConeMosaic.patternSampleSize(1);
        tBins = size(theEMPathsMeters,1);
        cMap = brewermap(tBins, 'Spectral');
        timeAxis = (0:1:(tBins-1))*theConeMosaic.integrationTime*1000;
        
        % Generate figure of the emPath
        figNo = 90;
        hFig = generateEMPathFig(figNo, theEMPathsMeters*1e6/theConeMosaic.micronsPerDegree*60, cMap, timeAxis);
        NicePlot.exportFigToPDF('componentFigs/EMpath.pdf', hFig, 300);

        % Generate figure of the cone mosaic with EMPath superimposed
        generateMosaicWithEMPathFig(theConeMosaic, theEMPathsMeters, tBins, cMap);
        
        % Generate figure of the photocurrent impulse response and
        % photocurrent noise
        tBinsIdx = size(photoCurrentsMean,4)+(-length(pCurrentImpulseResponses):-1);
        noiseTraces = bsxfun(@minus, photoCurrents(:,:,:,tBinsIdx), photoCurrentsMean(:,:,:,tBinsIdx));
        noiseTraces = reshape(noiseTraces, [size(photoCurrents,1)*size(photoCurrents,2)*size(photoCurrents,3) numel(tBinsIdx)]);
        % Only display 1024 traces
        noiseTraces = noiseTraces(1:1024,:);
        figNo = 91;
        hFig = generatePhotoCurrentFig(figNo, timeAxis, noiseTraces, pCurrentImpulseResponses);
        NicePlot.exportFigToPDF('componentFigs/Photocurrent.pdf', hFig, 300);
        
        % Visualize the first cone response mosaic instance at specific times together with the EMPath
        depictedTimes = 0:50:550;
        
        visualizedTrialConeExcitationResponseXYT = squeeze(coneExcitations(visualizedTrialIndex,:,:,:));
        
        videoFileName = 'response.mp4';
        videoOBJ = VideoWriter(videoFileName, 'MPEG-4'); % H264 format
        videoOBJ.FrameRate = 30;
        videoOBJ.Quality = 100;
        videoOBJ.open();
        
        for k = 1:numel(depictedTimes)
            [~,tBin] = min(abs((timeAxis-depictedTimes(k))));
            hFig = generateMosaicResponseWithEMPathFig(100+k,theConeMosaic, theEMPathsMeters, visualizedTrialConeExcitationResponseXYT, ...
                visualizedConePosMicrons, visualizedConeType, tBin, cMap, timeAxis(tBin));
            NicePlot.exportFigToPDF(sprintf('componentFigs/MosaicResponseWithEMPathComponent_%2.0fmsec.pdf',timeAxis(tBin)), hFig, 300);
            videoOBJ.writeVideo(getframe(hFig));
        end
        videoOBJ.close();
        fprintf('File saved in %s\n', videoFileName);
        
        
        % Collect all responses for the selected cones
        trialsNum = size(coneExcitations,1);
        visualizedConeExcitationResponses = zeros(3,trialsNum, tBins);
        visualizedPhotocurrentResponses = zeros(3,trialsNum, tBins);
        visualizedMeanConeExcitationResponse = zeros(3, tBins);
        visualizedMeanPhotocurrentResponse = zeros(3, tBins);
        
        for coneIdx = 1:3
            for trialIndex = 1:trialsNum
                tmp = squeeze(coneExcitations(trialIndex,:,:,:));
                tmp = reshape(tmp, [size(tmp,1)*size(tmp,2) size(tmp,3)]);
                visualizedConeExcitationResponses(coneIdx,trialIndex,:) = tmp(visualizedConeIndex(coneIdx),:);

                tmp = squeeze(photoCurrents(trialIndex,:,:,:));
                tmp = reshape(tmp, [size(tmp,1)*size(tmp,2) size(tmp,3)]);
                visualizedPhotocurrentResponses(coneIdx,trialIndex,:) = tmp(visualizedConeIndex(coneIdx),:);

                tmp = squeeze(coneExcitationsMean(1,:,:,:));
                tmp = reshape(tmp, [size(tmp,1)*size(tmp,2) size(tmp,3)]);
                visualizedMeanConeExcitationResponse(coneIdx,:) = tmp(visualizedConeIndex(coneIdx),:);

                tmp = squeeze(photoCurrentsMean(1,:,:,:));
                tmp = reshape(tmp, [size(tmp,1)*size(tmp,2) size(tmp,3)]);
                visualizedMeanPhotocurrentResponse(coneIdx,:) = tmp(visualizedConeIndex(coneIdx),:);
            end
        end
        
        % Visualize the responses of the selected cone
        figNo = 1000;
        hFig = generateSingleConeResponsePlots(figNo, timeAxis, visualizedTimeRange, ...
            visualizedMeanConeExcitationResponse, visualizedConeExcitationResponses, ...
            visualizedMeanPhotocurrentResponse, visualizedPhotocurrentResponses, visualizedConeType, cMap);
        NicePlot.exportFigToPDF('componentFigs/SingleConeResponses.pdf', hFig, 300);
           
    end       
end      
     
function hFig = generatePhotoCurrentFig(figNo, timeAxis, noiseTraces, pCurrentImpulseResponses)
    hFig = figure(figNo); clf
    set(hFig, 'Position', [10 10 426 420], 'Color', [1 1 1]);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 1, ...
       'heightMargin',  0.03, ...
       'widthMargin',    0.09, ...
       'leftMargin',     0.13, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.15, ...
       'topMargin',      0.02);
   
    ax = subplot('Position', subplotPosVectors(1,1).v);
    dt = timeAxis(2)-timeAxis(1);
    t = (0:(size(pCurrentImpulseResponses,1)-1))*dt;
    plot(t, pCurrentImpulseResponses(:,1), 'r-', 'Color', [1 0.2 0.4], 'LineWidth', 2.0); hold on
    plot(t, pCurrentImpulseResponses(:,2), 'g-', 'Color', [0.2 1.0 0.5], 'LineWidth', 2.0);
    plot(t, pCurrentImpulseResponses(:,3), 'b-', 'Color', [0.5 0.2 1.0], 'LineWidth', 2.0);
    
    plot(t, pCurrentImpulseResponses(:,1), 'r-', 'Color', [0 0 0], 'LineWidth', 3.0);
    plot(t, pCurrentImpulseResponses(:,2), 'g-', 'Color', [0 0 0], 'LineWidth', 3.0);
    plot(t, pCurrentImpulseResponses(:,3), 'b-', 'Color', [0 0 0], 'LineWidth', 3.0);
    
    plot(t, pCurrentImpulseResponses(:,1), 'r-', 'Color', [1 0.2 0.4], 'LineWidth', 2.0);
    plot(t, pCurrentImpulseResponses(:,2), 'g-', 'Color', [0.2 1.0 0.5], 'LineWidth', 2.0);
    plot(t, pCurrentImpulseResponses(:,3), 'b-', 'Color', [0.5 0.2 1.0], 'LineWidth', 2.0);
    legend({'L-cone', 'M-cone', 'S-cone'})
   

    grid on;
    set(gca, 'XLim', [t(1) t(end)]); 
    set(gca, 'XTick', [0:100:1000], 'YTick', [], 'YColor', 'none', 'XColor', 'none', 'XTickLabel', {});
    set(gca, 'FontSize', 20);
    box off
    
    ax = subplot('Position', subplotPosVectors(2,1).v);
    hPlot = plot(timeAxis(1:size(noiseTraces,2)),  noiseTraces, 'k-', 'LineWidth', 2);
    for k = 1:numel(hPlot)
       hPlot(k).Color(4) = 0.02;  % 5% transparent
    end
    
    grid on;
    set(gca, 'XLim', [t(1) t(end)]); 
    set(gca, 'XTick', [0:100:1000], 'YTick', [-10:5:10], 'YLim', [-9 9]);
    xlabel('time (msec)');
    ylabel('pAmps');
    set(gca, 'FontSize', 20);
    box off
    
end


function hFig = generateSingleConeResponsePlots(figNo, timeAxis, visualizedTimeRange, ...
    visualizedMeanConeExcitationResponse, visualizedConeExcitationResponses, ...
    visualizedMeanPhotocurrentResponse, visualizedPhotocurrentResponses, visualizedConeType, cMap)

    hFig = figure(figNo); clf
    set(hFig, 'Position', [10 10 780 1050], 'Color', [1 1 1]);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 3, ...
       'colsNum', 2, ...
       'heightMargin',  0.03, ...
       'widthMargin',    0.09, ...
       'leftMargin',     0.05, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.07, ...
       'topMargin',      0.02);
   
    XLim = visualizedTimeRange;
    XTick = 0:100:600;
    YTickConeExcitations = 0:10:100;
    YLimConeExcitations = [-2 85];
    YTickPhotocurrents = -80:10:0;
    YLimPhotocurrents = [-82 -8];
    
    for coneIdx = 1:numel(visualizedConeType)
        % Render cone excitation responses at the top
        ax = subplot('Position', subplotPosVectors(coneIdx,1).v);
        plotMeanAndSingleResponseInstances(ax, timeAxis, squeeze(visualizedMeanConeExcitationResponse(coneIdx,:)), ...
            squeeze(visualizedConeExcitationResponses(coneIdx,:,:)), ...
            XLim, YLimConeExcitations, XTick, YTickConeExcitations, 'cone excitations (R*/c/5 msec)', ...
            cMap, visualizedConeType(coneIdx), coneIdx == 3,coneIdx == 1);

        % Render pCurrent responses at the bottom
        ax = subplot('Position', subplotPosVectors(coneIdx,2).v);
        plotMeanAndSingleResponseInstances(ax, timeAxis, squeeze(visualizedMeanPhotocurrentResponse(coneIdx,:)), ...
            squeeze(visualizedPhotocurrentResponses(coneIdx,:,:)), ...
            XLim, YLimPhotocurrents, XTick, YTickPhotocurrents, 'photocurrent (pAmps)', ...
            cMap, visualizedConeType(coneIdx), coneIdx == 3, coneIdx == 1);
    end
    
    
end


function plotMeanAndSingleResponseInstances(ax, timeAxis, meanResponse, responseInstances, ...
            XLim, YLim, XTick, YTick, yAxisLabel, cMap, visualizedConeType, displayXLabel, displayTitle)
        
    tBins = numel(timeAxis);
    plot(timeAxis, responseInstances, 'k-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2); hold on
    plot(ax, timeAxis, meanResponse, 'k-', 'LineWidth', 6);
    % Overlay the trace in color
    for k = 1:tBins-1
       plot([timeAxis(k) timeAxis(k+1)], [meanResponse(k) meanResponse(k+1)], ...
          '-', 'Color', squeeze(cMap(k,:)), 'LineWidth', 4);
    end
    set(gca, 'XLim', XLim, 'XTick', XTick, 'YLim', YLim, 'YTick', YTick, 'FontSize', 20, 'LineWidth', 1.0);
    grid on; box on;
    switch(visualizedConeType-1)
        case 1
            coneType = 'L-cone';
            color = [1 0.5 0.4];
        case 2
            coneType = 'M-cone';
            color = [0.2 1 0.5];
        case 3
            coneType = 'S-cone';
            color = [0.5 0.2 1.0];
    end

        
    if (displayXLabel)
        xlabel('\it time (msec)');
    else
        set(gca, 'XTickLabel', {});
    end
    
    if (displayTitle)
        title(sprintf('%s', yAxisLabel));
    end
    if (strcmp(yAxisLabel, 'cone excitations (R*/c/5 msec)'))
        text(415, 75, coneType, 'Color', color, 'FontSize', 24, 'FontWeight', 'bold');
    end
    
    
end


function hFig = generateMosaicResponseWithEMPathFig(figNo, theConeMosaic, theEMPathsMeters, visualizedTrialConeExcitationResponseXYT, ...
    visualizedConePosMicrons, visualizedConeType, tBin, cMap, timeMsec)
    hFig = figure(figNo); clf
    set(hFig, 'Position', [10 10 426 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    
    theConeMosaic.renderActivationMap(ax, squeeze(visualizedTrialConeExcitationResponseXYT(:,:,tBin)),...
                'outlineConesAlongHorizontalMeridian', ~true, ...
                'mapType', 'modulated disks', ...
                'visualizedConeAperture', 'geometricArea');
    hold on;
    for coneIdx = 1:3
        switch(visualizedConeType(coneIdx)-1)
            case 1
                color = [1 0.5 0.4];
            case 2
                color = [0.2 1 0.5];
            case 3
                color = [0.5 0.2 1.0];
        end
        plot(ax,visualizedConePosMicrons(coneIdx,1)*1e-6, visualizedConePosMicrons(coneIdx,2)*1e-6, 'o', ...
            'Color', color, 'MarkerFaceColor', color, 'LineWidth', 1.0, 'MarkerSize', 8);
    end
    
    plotCrossHairs(ax, theEMPathsMeters, tBin, cMap);
    textPosMicrons = [-82 80];
    text(ax, textPosMicrons(1)*1e-6, textPosMicrons(2)*1e-6, sprintf('%2.0f msec', timeMsec), 'Color', 'w', 'FontSize', 26);
    hold off;
    set(ax,'XColor', 'none', 'YColor', 'none');
    set(ax,'XTickLabels', {});
    set(ax,'YTickLabels', {});
end

 
function hFig = generateEMPathFig(figNo, theEMPath, cMap, timeAxis)
    hFig = figure(figNo); clf
    set(hFig, 'Position', [10 10 426 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.11 0.16 0.87 0.82]);
    
    plot(ax, timeAxis, theEMPath(:, 1) , 'k-', 'LineWidth', 2); hold on;
    plot(ax, timeAxis, theEMPath(:, 2) , 'k--', 'LineWidth', 2); hold on;
    
    plot(ax, timeAxis, theEMPath(:, 1) , 'k-', 'LineWidth', 5); hold on;
    plot(ax, timeAxis, theEMPath(:, 2) , 'k:', 'LineWidth', 5); hold on;
    
    % Plot the trace in color
    for k = 1:(numel(timeAxis)-1)
        plot(ax, [timeAxis(k) timeAxis(k+1)], [theEMPath(k, 1) theEMPath(k+1, 1)], ...
            '-', 'Color', squeeze(cMap(k,:)), 'LineWidth', 2);
        plot(ax, [timeAxis(k) timeAxis(k+1)], [theEMPath(k, 2) theEMPath(k+1, 2)], ...
            '-', 'Color', squeeze(cMap(k,:)), 'LineWidth', 2);
    end
    legend({'x-pos', 'y-pos'}, 'Location', 'NorthWest');
    grid on; box on;
    set(gca, 'FontSize', 20, 'LineWidth', 1.0, 'YLim', [-5 17]);
    xlabel('\it time (msec)');
    ylabel('\it space (arc min)');
    axis 'square';
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
                'visualizedConeAperture', 'geometricArea', ...
                'apertureShape', 'disks', ...
                'backgroundColor', [1 1 1]);
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
        'apertureShape', 'disks', ...
        'backgroundColor', [1 1 1], ...
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

function theConeMosaic = generateMosaicFig(rootPath, mosaicFOV, theConeMosaic)

    if (isempty(theConeMosaic))
        load(fullfile(rootPath, sprintf('coneMosaic_%1.2fdegFOV.mat', mosaicFOV)), 'theConeMosaic');
    end
    
    hFig = figure(4); clf
    set(hFig, 'Position', [10 10 426 420], 'Color', [1 1 1]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    
    theConeMosaic.visualizeGrid('axesHandle', ax, ...
        'visualizedConeAperture', 'geometricArea', ...
        'apertureShape', 'disks', ...
        'backgroundColor', [1 1 1]);
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


