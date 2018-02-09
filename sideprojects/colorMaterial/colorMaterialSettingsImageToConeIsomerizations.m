function colorMaterialSettingsImageToConeIsomerizations

    close all
    
    [localDir,~] = fileparts(which(fullfile(mfilename)));
    
    horizontalFOV = 15;
    maxGridAdjustmentIterations = 50;
    
    regenerateScenes = ~true;
    regenerateConeMosaic = ~true;
    generateOI = ~true;
    generateResponses = true;
    
    visualizeScenes = ~true;
    visualizeMosaic = ~true;
    visualizeSceneOIandMosaic = ~true;

    
    % Mat filename containing previously computed scenes
    scenesMatFileName = sprintf('scenes.mat');
    
    % Mat filename containing previously computed cone mosaic
    coneMosaicMatFileName = sprintf('coneMosaic_%2.2fdegFOV.mat', horizontalFOV);

    % Mat filename containing previously computed oi
    oiMatFileName = sprintf('oi_%2.2fdegFOV.mat', horizontalFOV);
    
    % Mat filename containing previously computed responses
    responsesMatFileName = sprintf('responses_%2.2fdegFOV.mat', horizontalFOV);
    
    
    if (generateResponses)
        % Generate scenes
        if (~regenerateScenes)
            load(scenesMatFileName, 'scenesList', 'noBorderSceneRows', 'noBorderSceneCols');   
        else
            % Stuff we need to generate the scenes
            tbUse({'BrainardLabBase', 'isetbio'})
            cd(localDir);

            % Generate display object corresponding to calFile used in ColorMaterial experiments
            [display,cal] = generateISETbioDisplayForColorMaterialExperiment('EyeTrackerLCD');
            [scenesList, noBorderSceneRows, noBorderSceneCols] = generateISETbioScenes(display);
            save(scenesMatFileName, 'scenesList', 'noBorderSceneRows', 'noBorderSceneCols', '-v7.3');
        end

        % Display the scenes
        if (visualizeScenes)
            displayScenes(scenesList, noBorderSceneRows, noBorderSceneCols);
        end
    end
    
    % Now we need some stuff from IBIOColorDetect
    %tbUseProject('IBIOColorDetect')
    cd(localDir);
    
    % Generate the cone mosaic
    if (regenerateConeMosaic)
        integrationTime = 100/1000;
        theConeMosaic = coneMosaicGenerate(horizontalFOV,integrationTime, maxGridAdjustmentIterations);
        theConeMosaic.visualizeGrid(...
            'apertureShape', 'disks', ...
            'visualizedConeAperture', 'lightCollectingArea', ...
            'labelConeTypes', true,...
            'generateNewFigure', true);
        save(coneMosaicMatFileName, 'theConeMosaic', '-v7.3');
    else
        load(coneMosaicMatFileName, 'theConeMosaic');
    end
    
    if (visualizeMosaic)
        theConeMosaic.visualizeGrid(); % 'visualizedConeAperture', 'geometricArea');
    end
      
    if (generateResponses)
        if (generateOI)
            % Generate human optics
            oiParams = struct(...
                'opticsModel', 'WvfHuman', ...
                'wavefrontSpatialSamples', 301, ...
                'pupilDiamMm', 3.0, ...
                'umPerDegree', 300);
            theOI = oiWithCustomOptics(oiParams.opticsModel, oiParams.wavefrontSpatialSamples, oiParams.pupilDiamMm, oiParams.umPerDegree);

            % Set the FOV
            theOI = oiSet(theOI,'h fov',horizontalFOV);

            % Set the fNumber
            focalLength = oiGet(theOI,'distance');
            desiredFNumber = focalLength/(oiParams.pupilDiamMm/1000);
            theOI  = oiSet(theOI ,'optics fnumber',desiredFNumber);

            % Save the OI
            save(oiMatFileName, 'theOI', '-v7.3');
        else
            load(oiMatFileName, 'theOI');
        end
    
        % Container to hold all sceneData
        allScenesData = containers.Map();

        % Loop through images
        for sceneIndex = 1:numel(scenesList)
            % Grab the scene from the list of scenes
            theScene = scenesList{sceneIndex};
            theSceneName = sceneGet(theScene, 'name');

            % Compute the OIimage for this scene
            fprintf('Computing optical image for scene %d\n', sceneIndex);
            theOI = oiCompute(theOI, theScene);

            % Visualize scene and optical image
            if (visualizeSceneOIandMosaic)
                visualizeStimulusAndConeMosaic(theConeMosaic, theOI, theScene, [])
            end

            % Compute the isomerizations for this optical image
            fprintf('Computing mosaic response for scene %d\n', sceneIndex);
            theFullMosaicIsomerizations = theConeMosaic.compute(theOI, 'currentFlag',false);

            % Isomerization rates
            theFullMosaicIsomerizations = theFullMosaicIsomerizations / theConeMosaic.integrationTime;
            maxIsomerizationRate = max(theFullMosaicIsomerizations(:));
               
            % Containers to hold various response components
            theConeIsomerizations = containers.Map();
            theDemosaicedIsomerizationsMaps = containers.Map();
            demosaicedIsomerizationsMapSupportDegs = containers.Map();
            xConeLocsDegs = containers.Map();
            yConeLocsDegs = containers.Map();
            
            % Compute responses for individual cone types
            coneTypeNames = {'L-cones', 'M-cones', 'S-cones', 'All-cones'};
            for coneType = 1:numel(coneTypeNames)
                coneName = coneTypeNames{coneType};
                demosaicingSampleSpacingMicrons = 1;
                [theDemosaicedIsomerizationsMaps(coneName), demosaicedIsomerizationsMapSupportDegs(coneName), ...
                    theConeIsomerizations(coneName), xConeLocsDegs(coneName), yConeLocsDegs(coneName)] = ...
                    theConeMosaic.demosaicConeTypeActivationFromFullActivation(coneName, ...
                         theFullMosaicIsomerizations, demosaicingSampleSpacingMicrons);
            end

            allSceneData = struct(...
                'theConeIsomerizations', theConeIsomerizations, ... 
                'xConeLocsDegs', xConeLocsDegs, ...                                 
                'yConeLocsDegs', yConeLocsDegs, ...                                
                'theDemosaicedIsomerizationsMaps', theDemosaicedIsomerizationsMaps, ... 
                'theDemosaicedIsomerizationsMapSupport', demosaicedIsomerizationsMapSupportDegs...
                );
            allSceneData.maxIsomerizationRate = maxIsomerizationRate;
            
            allScenesData(theSceneName) = allSceneData;
        end % sceneIndex
    
        fprintf('Finished with everything.\nExporting responses to %s....', responsesMatFileName);
        save(responsesMatFileName, 'allScenesData', '-v7.3');
        fprintf('Finished exporting responses.\n');
    end
end
 
function cm = coneMosaicGenerate(fovDegs, integrationTime, maxGridAdjustmentIterations)

    quality.resamplingFactor = 2;
    quality.tolerance1 = 1.0;           % larger than default tolerances to speed-up computation. For production work, either do not set, or set to equal or lower than 0.01 
    quality.tolerance2 = 0.5;           % larger than default tolerances to speed-up computation, For production work, either do not set, or set to equal or lower than 0.001 
    quality.marginF = [];    
        
    % Instantiate a hex mosaic
    cm = coneMosaicHex(quality.resamplingFactor, ...
        'fovDegs', fovDegs, ...
        'eccBasedConeDensity', true, ... 
        'maxGridAdjustmentIterations', maxGridAdjustmentIterations, ...
        'latticeAdjustmentPositionalToleranceF', quality.tolerance1, ...         
        'latticeAdjustmentDelaunayToleranceF', quality.tolerance2, ...     
        'marginF', quality.marginF ...
    );

    % Set the mosaic's integration time
    cm.integrationTime = integrationTime;
    cm.noiseFlag = 'none';
end

function displayScenes(scenesList, noBorderSceneRows, noBorderSceneCols)
    
    rows = 3; cols = 5;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rows, ...
           'colsNum', cols, ...
           'heightMargin',   0.04, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.04, ...
           'topMargin',      0.03);
       
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1500 940]);
    for k = 1:numel(scenesList)
        sceneName = sceneGet(scenesList{k}, 'name');
        sceneLuminance = sceneGet(scenesList{k}, 'luminance');
        sceneMeanLuminance = mean(mean(sceneLuminance(noBorderSceneRows,noBorderSceneCols)));
        sceneHorizontalFOV = sceneGet(scenesList{k}, 'wangular');
        sceneVerticalFOV = sceneGet(scenesList{k}, 'hangular');
        sceneRGB = sceneGet(scenesList{k}, 'rgb image');
        angRes = sceneGet(scenesList{k}, 'angular resolution');
        
        xspace = (1:size(sceneRGB,2))*angRes(1);
        xspace = xspace - mean(xspace);
        yspace = (1:size(sceneRGB,2))*angRes(2);
        yspace = yspace - mean(yspace);
        
        row = floor((k-1)/cols)+1;
        col = mod(k-1, cols)+1;
        subplot('Position', subplotPosVectors(row,col).v);
        image(xspace, yspace, sceneRGB);
        title(sprintf('%s\nmean lum=%2.1f cd/m2, FOV=%2.1f x %2.1f deg', sceneName, sceneMeanLuminance, sceneHorizontalFOV, sceneVerticalFOV));
        axis 'xy'
        axis 'image';
 
        set(gca, 'FontSize', 12);
        if (row == rows) && (col == 1)
            xlabel('space (degs)');
            ylabel('space (degs)');
        else
            set(gca, 'XTickLabel', {}, 'YTickLabel', {});
        end
        
    end
    NicePlot.exportFigToPNG('scenes.png', hFig,300);
end

function [sceneList, noBorderSceneRows, noBorderSceneCols] = generateISETbioScenes(displayOBJ)

    %% Specify images to load
    RGBsettingsImageDir = 'ColorMatStimuliRGB';
    imageFileNames = {...
        'C1M4' ...
        'C2M4' ...
        'C3M4' ...
        'C4M1' ...
        'C4M2' ...
        'C4M3' ...
        'C4M4' ...
        'C4M5' ...
        'C4M6' ...
        'C4M7' ...
        'C5M4' ...
        'C6M4' ...
        'C7M4' ...
        };
    
    for imageIndex = 1:numel(imageFileNames)
        
        %% Load RGB settings image
        load(fullfile(RGBsettingsImageDir, sprintf('%s-RGB.mat',imageFileNames{imageIndex})), 'sensorImageRGB');
      
        %% Add border
        [rows, cols, ~] = size(sensorImageRGB);
        borderHeight = floor(cols-rows)/2;
        borderWidth = 0;
        
        extraBorder = 0;
        borderWidth = borderWidth + extraBorder;
        borderHeight = borderHeight + extraBorder;
        
        if (borderWidth > 0) || (borderHeight > 0)
            [rows, cols, ~] = size(sensorImageRGB);
            tmp = sensorImageRGB;

            if (imageIndex == 1)
                borderRGB = 0*tmp(1,1,:);
            end

            sensorImageRGB = bsxfun(@plus, zeros(rows+2*borderHeight, cols+2*borderWidth, 3), borderRGB);
            noBorderSceneRows = borderHeight+(1:rows);
            noBorderSceneCols = borderWidth+(1:cols);
            for channel = 1:3
                theSlice = squeeze(tmp(:,:,channel));
                sensorImageRGB(noBorderSceneRows, noBorderSceneCols, channel) = flipud(theSlice);
            end
        end
        
        %% Convert to ISETbio scene
        meanLuminance = [];
        scene = sceneFromFile(sensorImageRGB, 'rgb', meanLuminance, displayOBJ);
        
        % Add a name to the scene
        scene = sceneSet(scene, 'name', imageFileNames{imageIndex});
        
        % Add to list of scenes
        sceneList{imageIndex} = scene;
    end
end

function [displayOBJ, cal] = generateISETbioDisplayForColorMaterialExperiment(calFileName)
    %% Load calibration file
    cal = LoadCalFile(calFileName,11);
    
    %% Subsample spectral axis to a 10 nm wavelength spacing
    wavelengths = SToWls(cal.describe.S);
    wavelengths = wavelengths(1):10:wavelengths(end);
    
    %% Generate ISETbio display object corresponding to calFile
    extraData = ptb.ExtraCalData;
    extraData.distance = 0.7;
    extraData.subSamplingSvector = WlsToS(wavelengths');
    saveDisplayObject = false;
    displayOBJ = ptb.GenerateIsetbioDisplayObjectFromPTBCalStruct(calFileName, cal, extraData, saveDisplayObject);
    % The effective pixel size models the fact that images were scaled
    % down, not presented at a 1:1 image pixel:screen pixel. According to
    % Ana Radonjic, the effective pixels size was 0.01875
    effectivePixelSizeCentiMeters = 0.01875;
    % Adjust DPI to account for this effective pixel size
    actualPixelSizeCentimeters = 100/displayGet(displayOBJ, 'dots per meter');
    actualDPI = displayGet(displayOBJ, 'dpi');
    effectiveDPI = actualDPI * actualPixelSizeCentimeters/effectivePixelSizeCentiMeters;
    displayOBJ = displaySet(displayOBJ, 'dpi', effectiveDPI);
    
    % Also return the cal in XYZ color space (for debugging purposes only)
    load T_xyz1931
    T_xyz1931 = 683*T_xyz1931;
    cal = SetSensorColorSpace(cal,T_xyz1931,S_xyz1931);    
end
