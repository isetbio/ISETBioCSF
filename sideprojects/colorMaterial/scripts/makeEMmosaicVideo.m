function makeEMmosaicVideo

    close all
    
    [localDir,~] = fileparts(which(fullfile(mfilename)));
    
    horizontalFOV = 2.0;
    visualizedFOVdegs = 0.4;
    
    showSceneOnBackground = false;
    showMovingMosaicOnSeparateSubFig = false;
    showCircle = false;
    showCross = true;
    showEMpathOnTopOfConeMosaic = true;
        
    % Mat filename containing previously computed cone mosaic
    coneMosaicMatFileName = sprintf('../coneMosaic_%2.2fdegFOV.mat', horizontalFOV);
    load(coneMosaicMatFileName, 'theConeMosaic');
    
    if (showSceneOnBackground)
        % Load the scenes
        scenesMatFileName = sprintf('../scenes.mat');
        load(scenesMatFileName, 'scenesList', 'noBorderSceneRows', 'noBorderSceneCols');   
        whichScene = 'C4M4';
        sceneFOV = [];
        whichSceneIndex = [];
        for k = 1:numel(scenesList)
                if (strfind(whichScene,scenesList{k}.name))
                    whichSceneIndex = k;
                end
        end
        scenesList = {scenesList{whichSceneIndex}};
        sceneRGB = sceneGet(scenesList{1}, 'rgb image');
        angRes = sceneGet(scenesList{1}, 'angular resolution');
        xspace = (1:size(sceneRGB,2))*angRes(1);
        xspaceDegs = xspace - mean(xspace);
        yspace = (1:size(sceneRGB,2))*angRes(2);
        yspaceDegs = yspace - mean(yspace);
        
        % Scale the scene as if it were 10 times further that it actually was
        viewingDistanceScalar = 1/10;
        xspaceMeters = viewingDistanceScalar *xspaceDegs * theConeMosaic.micronsPerDegree * 1e-6;
        yspaceMeters = viewingDistanceScalar *yspaceDegs * theConeMosaic.micronsPerDegree * 1e-6;
    end
    
    
    theConeMosaic.integrationTime = 0.2/1000;
    fixEMobj = fixationalEM();
    fixEMobj.microSaccadeMeanIntervalSeconds = 0.250;
    
    eyeMovementsPerTrial = 1.001/theConeMosaic.integrationTime;
    nTrials = 1;
    
    fixEMobj.computeForConeMosaic(theConeMosaic, eyeMovementsPerTrial, ...
        'nTrials', nTrials, ...
        'rSeed', 857, ...
        'useParfor', true);
    
    crossLengthDegs = 0.05;
    crossLengthMeters = crossLengthDegs * theConeMosaic.micronsPerDegree * 1e-6;

   
    emRadiusDegs = 0.22;
    emRadiusMeters = emRadiusDegs * theConeMosaic.micronsPerDegree * 1e-6;
    emRegion.x = emRadiusMeters * cosd(0:5:360);
    emRegion.y = emRadiusMeters * sind(0:5:360);
        
    
    visualizedSpaceMeters = visualizedFOVdegs * theConeMosaic.micronsPerDegree * 1e-6;

    ticksDegs = -0.5:0.1:0.5;
    ticks = ticksDegs * theConeMosaic.micronsPerDegree * 1e-6;
    tickLabels = sprintf('%2.1f\n', ticksDegs);
        
    videoOBJ = VideoWriter('emMosaic.mp4', 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 60;
    videoOBJ.Quality = 100;
    videoOBJ.open();
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', 1+showMovingMosaicOnSeparateSubFig, ...
           'heightMargin',  0.08, ...
           'widthMargin',    0.05, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.02, ...
           'bottomMargin',   0.06, ...
           'topMargin',      0.02);
    
    if (showMovingMosaicOnSeparateSubFig)
        figSize = [1100 580];
    elseif (showEMpathOnTopOfConeMosaic)
        figSize = [740 770];
    else
        figSize = [510 580];
    end
    
    if (showEMpathOnTopOfConeMosaic)
        coneLocsMeters = theConeMosaic.coneLocsHexGrid;
        idx = sqrt(sum(coneLocsMeters .^2, 2)) < visualizedSpaceMeters/2;
        coneLocsMeters  = coneLocsMeters (idx,:);
        coneApertureRadiusMeters = diameterForCircularApertureFromWidthForSquareAperture(theConeMosaic.pigment.width)/2;
        coneOutline.x = coneApertureRadiusMeters * cosd(0:30:360);
        coneOutline.y = coneApertureRadiusMeters * sind(0:30:360);
    end
    
    for iTrial = 1:nTrials
        
        hFig = figure(1); clf
        set(hFig, 'Position', [10 10 figSize(1) figSize(2)], 'Color', [1 1 1]);
    
        if (showMovingMosaicOnSeparateSubFig)
            axMosaic = subplot('Position', subplotPosVectors(1,1).v);
            axEMpath = subplot('Position', subplotPosVectors(1,2).v);
        else
            axEMpath = subplot('Position', subplotPosVectors(1,1).v);
        end
        
        
        
        if (showMovingMosaicOnSeparateSubFig)
            theConeMosaic.visualizeGrid('axeshandle', axMosaic, 'backgroundColor', [.75 .75 .75], ...
                'foregroundColor', [0 0 0], ...
                'visualizedConeAperture', 'geometricArea', ...
                'labelConeTypes', true,...
                'apertureshape', 'disks');

            hold(axMosaic, 'on');
            plot(axMosaic,[0 0], crossLengthMeters*[-1 1], 'k-', 'LineWidth', 2.0);
            plot(axMosaic, crossLengthMeters*[-1 1], [0 0], 'k-', 'LineWidth', 2.0);
            hold(axMosaic, 'off');
            box(axMosaic, 'on');
        end
        
        % emPath for this trial
        theEMpathMeters = squeeze(fixEMobj.emPosMicrons(iTrial,:,:)) * 1e-6;
    
        for k = 1:size(theEMpathMeters,1)
            if (showMovingMosaicOnSeparateSubFig)
                xo = theEMpathMeters(k,1);
                yo = theEMpathMeters(k,2);
                set(axMosaic, 'XLim', xo + visualizedSpaceMeters*[-0.5 0.5]*1.05, 'YLim', yo + visualizedSpaceMeters*[-0.5 0.5]*1.05, 'XTickLabel', {}, 'YTickLabel', {});
                set(axMosaic, 'FontSize', 16, 'LineWidth', 1.0);
                title(axMosaic,sprintf('%2.1f msec', 1000*fixEMobj.timeAxis(k)));
            end
            
            if (showSceneOnBackground)
                image(axEMpath, xspaceMeters, yspaceMeters, sceneRGB);
            end
            if (showCircle)
                plot(axEMpath, emRegion.x, emRegion.y, '--', 'Color', [0 0 0], 'LineWidth', 1.0);
            end
            
            if (showEMpathOnTopOfConeMosaic)
                for iCone = 1:size(coneLocsMeters,1)
                    plot(coneLocsMeters(iCone,1)+coneOutline.x, coneLocsMeters(iCone,2)+coneOutline.y, 'k-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);
                    if (iCone == 1)
                        hold(axEMpath, 'on')
                    end
                end
            end
            
            hold(axEMpath, 'on')
            plot(axEMpath, theEMpathMeters(1:k,1),   theEMpathMeters(1:k,2), 'r-', 'Color', [1 0.5 0.5], 'LineWidth', 3.0);
            plot(axEMpath, theEMpathMeters(1:k,1),   theEMpathMeters(1:k,2), 'r-', 'LineWidth', 1.5);
            if (showCross)
                plot(axEMpath,theEMpathMeters(k,1)+[0 0],theEMpathMeters(k,2)+crossLengthMeters*[-1 1], 'k-', 'LineWidth', 1.5);
                plot(axEMpath,theEMpathMeters(k,1)+crossLengthMeters*[-1 1],  theEMpathMeters(k,2)+[0 0],  'k-', 'LineWidth', 1.5);
            end
            axis(axEMpath, 'xy'); axis(axEMpath, 'image');
            plot(axEMpath,[0 0], 1e-6*[-300 300], '-', 'Color', [0.3 0.3 0.3]);
            plot(axEMpath,1e-6*[-300 300],[0 0],  '-', 'Color', [0.3 0.3 0.3]);
            
            axis 'square'
            hold(axEMpath, 'off');
            set(axEMpath, 'XLim', visualizedSpaceMeters*[-0.5 0.5]*1.05, 'YLim', visualizedSpaceMeters*[-0.5 0.5]*1.05, ...
                'XTick', ticks, 'YTick', ticks, 'XTickLabels', tickLabels, 'YTickLabels', {}, ...
                'FontSize', 16);
            grid(axEMpath, 'on'); box(axEMpath, 'on');
            set(axEMpath, 'LineWidth', 1.0);
            if (nTrials>1)
                title(axEMpath,sprintf('%2.1f msec (trial #%d)',1000*fixEMobj.timeAxis(k), iTrial));
            else
                title(axEMpath,sprintf('%2.1f msec',1000*fixEMobj.timeAxis(k)));
            end
            xlabel(axEMpath, 'space (degs)');
            drawnow;
            videoOBJ.writeVideo(getframe(hFig));
        end
    end
    
     videoOBJ.close();
    
        
end

