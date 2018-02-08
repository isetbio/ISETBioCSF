function settingsImageToConeIsomerizations

    [localDir,~] = fileparts(which(fullfile(mfilename)));
    
    % Stuff we need to generate the scenes
    tbUse({'BrainardLabBase', 'isetbio'})
    cd(localDir);
    
    desiredHorizontalFOV = 2;
    maxGridAdjustmentIterations = 100;
    
    regenerateScenes = ~true;
    regenerateConeMosaic = ~true;
    generateOI = true;
    
    
    visualizeScenes = ~true;
    visualizeMosaic = ~true;
    visualizeSceneOIandMosaic = ~true;
    visualizeResponses = true;
    
    % Mat filename containing previously computed scenes
    scenesMatFileName = sprintf('scenes.mat');
    
    % Mat filename containing previously computed cone mosaic
    coneMosaicMatFileName = sprintf('coneMosaic_%2.2fdegFOV.mat', desiredHorizontalFOV);

    % Mat filename containing previously computed oi
    oiMatFileName = sprintf('oi_%2.2fdegFOV.mat', desiredHorizontalFOV);
    
    % Mat filename containing the computed responses
    responsesMatFileName = sprintf('responses_%2.2fdegFOV.mat', desiredHorizontalFOV);
    
    %% Load scenes
    if (~regenerateScenes)
        load(scenesMatFileName, 'scenesList', 'noBorderSceneRows', 'noBorderSceneCols');   
    else
        %% Generate display object corresponding to calFile used in ColorMaterial experiments
        [display,cal] = generateISETbioDisplayForColorMaterialExperiment('EyeTrackerLCD');
        cal = [];
        [scenesList, noBorderSceneRows, noBorderSceneCols] = generateISETbioScenes(display, cal);
        save(scenesMatFileName, 'scenesList', 'noBorderSceneRows', 'noBorderSceneCols', '-v7.3');
    end
    
    %% Display the scenes
    if (visualizeScenes)
        displayScenes(scenesList, noBorderSceneRows, noBorderSceneCols);
    end
    
    
    % Now we need some stuff from IBIOColorDetect
    tbUseProject('IBIOColorDetect')
    cd(localDir);
    
    %% Generate the cone mosaic
    if (regenerateConeMosaic)
        integrationTime = 100/1000;
        theConeMosaic = coneMosaicGenerate(desiredHorizontalFOV,integrationTime, maxGridAdjustmentIterations);
        theConeMosaic.visualizeGrid(...
            'apertureShape', 'disks', ...
            'visualizedConeAperture', 'lightCollectingArea', ...
            'labelConeTypes', true,...
            'generateNewFigure', true);
        save(coneMosaicMatFileName, 'theConeMosaic', '-v7.3');
    else
        load(coneMosaicMatFileName, 'theConeMosaic');
    end
    
    
    if (generateOI)
        %% Generate human optics
        oiParams = struct(...
            'opticsModel', 'WvfHuman', ...
            'wavefrontSpatialSamples', 301, ...
            'pupilDiamMm', 3.0, ...
            'umPerDegree', 300);
        theOI = oiWithCustomOptics(oiParams.opticsModel, oiParams.wavefrontSpatialSamples, oiParams.pupilDiamMm, oiParams.umPerDegree);

        %% Set the FOV
        theOI = oiSet(theOI,'h fov',desiredHorizontalFOV);

        %% Set the fNumber
        focalLength = oiGet(theOI,'distance');
        desiredFNumber = focalLength/(oiParams.pupilDiamMm/1000);
        theOI  = oiSet(theOI ,'optics fnumber',desiredFNumber);
        
        %% Save the OI
        save(oiMatFileName, 'theOI', '-v7.3');
    else
        load(oiMatFileName, 'theOI');
    end
    
    if (visualizeMosaic)
        theConeMosaic.visualizeGrid(); % 'visualizedConeAperture', 'geometricArea');
    end
    
    allSceneData = containers.Map();
    
    %% Loop through images
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
        
        % Compute the isomerization rate for this optical image
        fprintf('Computing mosaic response for scene %d\n', sceneIndex);
        theIsomerizations = theConeMosaic.compute(theOI, 'currentFlag',false) / theConeMosaic.integrationTime;
        
        % Extract the full-pattern isomerizations
        theIsomerizations = reshape(theIsomerizations, [prod(size(theIsomerizations)) 1]);
        
        % Demosaic response
        demosaicingSampleSpacingMicrons = 1;
        demosaicRangeMicrons = max(theConeMosaic.fov)*theConeMosaic.micronsPerDegree;
        
        coneTypeNames = {'L-cones', 'M-cones', 'S-cones', 'All-cones'};
        maxIsomerizationRate = 0;
        theConeIsomerizations = containers.Map();
        theDemosaicedIsomerizationsMaps = containers.Map();
        xConeLocsDegs = containers.Map();
        yConeLocsDegs = containers.Map();
        
        for coneType = 1:numel(coneTypeNames)
            % Compute demosaiced-map
            coneName = coneTypeNames{coneType};
            [theDemosaicedIsomerizationsMap, support, theConeIsomerizations(coneName), xConeLocsDegs(coneName), yConeLocsDegs(coneName)] = ...
                demosaicResponseFromFullPatternResponse(coneName, ...
                     theConeMosaic, theIsomerizations, demosaicingSampleSpacingMicrons, demosaicRangeMicrons);
            % Isomerization rate
            theDemosaicedIsomerizationsMap = theDemosaicedIsomerizationsMap / theConeMosaic.integrationTime;
            theConeIsomerizations(coneName) = theConeIsomerizations(coneName) / theConeMosaic.integrationTime;
            
            m = max(theDemosaicedIsomerizationsMap(:));
            if (m > maxIsomerizationRate)
                maxIsomerizationRate = m;
            end
            theDemosaicedIsomerizationsMaps(coneName) = theDemosaicedIsomerizationsMap;
        end
        
        allSceneData(theSceneName) = struct(...
            'theScene', theScene, ...
            'theOI', theOI, ...
            'coneTypeNames', coneTypeNames, ...
            'theConeIsomerizations', theConeIsomerizations, ...                 % for 4 cone types
            'xConeLocsDegs', xConeLocsDegs, ...                                 % for 4 cone types
            'yConeLocsDegs', yConeLocsDegs, ...                                 % for 4 cone types
            'theDemosaicedIsomerizationsMaps', theDemosaicedIsomerizationsMaps, ...  % for 4 cone types
            'theDemosaicedIsomerizationsMapSupport', 'support' ...
            );
        
        if (visualizeResponses)
            hFig = figure(25); clf;
            set(hFig, 'Position', [10 10 1450 780]);
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', 2, ...
               'colsNum', numel(coneTypeNames), ...
               'heightMargin',   0.04, ...
               'widthMargin',    0.02, ...
               'leftMargin',     0.01, ...
               'rightMargin',    0.006, ...
               'bottomMargin',   0.02, ...
               'topMargin',      0.02);

            for coneType = 1:numel(coneTypeNames)
                coneName = coneTypeNames{coneType};
                theDemosaicedIsomerizationsMap = theDemosaicedIsomerizationsMaps(coneName);
                
                subplot('Position', subplotPosVectors(1,coneType).v);
                imagesc(support,support,theDemosaicedIsomerizationsMap/maxIsomerizationRate);
                set(gca, 'CLim', [0 1], 'XLim', [support(1) support(end)], 'YLim', [support(1) support(end)]);
                axis 'image'; axis 'xy';
                title(sprintf('%s (max: %2.0f R*/sec', coneTypeNames{coneType}, max(theDemosaicedIsomerizationsMap(:))));

                subplot('Position', subplotPosVectors(2,coneType).v);
                imagesc(support,support,theDemosaicedIsomerizationsMap/max(theDemosaicedIsomerizationsMap(:)));
                set(gca, 'CLim', [0 1], 'XLim', [support(1) support(end)], 'YLim', [support(1) support(end)]);
                axis 'image'; axis 'xy';
                title(sprintf('%s (max: %2.0f R*/sec', coneTypeNames{coneType}, max(theDemosaicedIsomerizationsMap(:))));
            end
            colormap(bone(1024));
        end
        
    end % sceneIndex
    
   save(responsesMatFileName, 'allSceneData', '-v7.3');
   fprintf('Finished with everything. Responses saved in %s.\n', responsesMatFileName);
    
end

function [demosaicedActivationMap, support, theConeTypeResponse, xConeLocsDegs, yConeLocsDegs] = demosaicResponseFromFullPatternResponse(coneType,...
    theConeMosaic, theFullPatternResponse, demosaicingSampleSpacingMicrons, demosaicRangeMicrons)
        
    xPatternSupport = squeeze(theConeMosaic.patternSupport(:,:,1));
    yPatternSupport = squeeze(theConeMosaic.patternSupport(:,:,2));
    
    switch coneType
        case 'L-cones'
            targetConesIndices = find(theConeMosaic.pattern==2);
        case 'M-cones'
            targetConesIndices = find(theConeMosaic.pattern==3);
        case 'S-cones'
            targetConesIndices = find(theConeMosaic.pattern==4);
        case 'All-cones'
            targetConesIndices = find(theConeMosaic.pattern > 1);
        otherwise
            error('Invalid cone type: ''%s''.', coneType);
    end
    
    theConeTypeResponse = theFullPatternResponse(targetConesIndices);
    xPatternSupport = xPatternSupport(targetConesIndices);
    yPatternSupport = yPatternSupport(targetConesIndices);
    
    xConeLocsDegs = xPatternSupport*1e6/theConeMosaic.micronsPerDegree;
    yConeLocsDegs = yPatternSupport*1e6/theConeMosaic.micronsPerDegree;
        
    degIncrement = demosaicingSampleSpacingMicrons/theConeMosaic.micronsPerDegree;
    demosaicRangeDegs = demosaicRangeMicrons / theConeMosaic.micronsPerDegree;
    support = -demosaicRangeDegs/2:degIncrement:demosaicRangeDegs/2;
    [X,Y] = meshgrid(support,support);
        
    % Compute demosaic response
    method = 'nearest'; %'natural'; %'nearest';
    demosaicedActivationMap = griddata(xConeLocsDegs, yConeLocsDegs, theConeTypeResponse(:), X,Y, method);
end

    
function cm = coneMosaicGenerate(fovDegs, integrationTime, maxGridAdjustmentIterations)

    quality.resamplingFactor = 2;
    quality.tolerance1 = 1.0;                           % larger than default tolerances to speed-up computation. For production work, either do not set, or set to equal or lower than 0.01 
    quality.tolerance2 = 0.5;                          % larger than default tolerances to speed-up computation, For production work, either do not set, or set to equal or lower than 0.001 
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

function [sceneList, noBorderSceneRows, noBorderSceneCols] = generateISETbioScenes(displayOBJ, cal)

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
    
    debug = false;

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
        
        if (debug)&&(~isempty(cal))
            % Only used to debug conversion to ISETbio
            [sensorImageRGBCalFormat, nCols, mRows] = ImageToCalFormat(sensorImageRGB);
            sensorImageXYZCalFormat = SettingsToSensor(cal,sensorImageRGBCalFormat);
            sensorImageXYZ = CalFormatToImage(sensorImageXYZCalFormat, nCols, mRows);
        end
        
        %% Convert to ISETbio scene
        meanLuminance = [];
        scene = sceneFromFile(sensorImageRGB, 'rgb', meanLuminance, displayOBJ);
        
%         if (imageIndex == 1)
%             if (debug)
%                 meanLumBefore = sceneGet(scene, 'mean lum')
%                 vertFOVbefore = sceneGet(scene, 'hangular')
%                 horizFOVbefore = sceneGet(scene, 'wangular')
%             end
%             
%             % Correct the dpi for the display so as to achieve the desired
%             % image FOV. We are doing this because the images were not
%             % presented at a 1:1 image:pixel ratio (they were shunk)
%             displayDotsPerCentiMeter = 100/displayGet(displayOBJ, 'dots per meter')
%             displayOBJ = displaySet(displayOBJ, 'dots per meter', 0.01875);
%             displaySizeMeters = displayGet(displayOBJ,'size')
%             dotSizeCentiMeters = 100/displayDotsPerMeter;
%             pixelSh
%             pause
%             
%             horizFOVforUnityPixelRatio = sceneGet(scene, 'wangular');
%             pixelShrinkingFactor = desiredHorizontalFOV/horizFOVforUnityPixelRatio;
%             DPI = displayGet(displayOBJ, 'dpi');
%             adjustedDPI = DPI/pixelShrinkingFactor;
%             displayOBJ = displaySet(displayOBJ, 'dpi', adjustedDPI);
%             fprintf('Adjusted display DPI from %2.3f to %2.3f to obtain desired stimulus FOV\n',DPI, adjustedDPI)
%             
%             % Get the scene based on the adjusted display
%             scene = sceneFromFile(sensorImageRGB, 'rgb', meanLuminance, displayOBJ);
%             
%             if (debug)
%                 meanLumAfter = sceneGet(scene, 'mean lum')
%                 vertFOVAfter = sceneGet(scene, 'hangular')
%                 horizFOVAfter = sceneGet(scene, 'wangular')
%             end
%         end % (imageIndex == 1)
        
        % Add a name tp the scene
        scene = sceneSet(scene, 'name', imageFileNames{imageIndex});
        
        % Add to list of scenes
        sceneList{imageIndex} = scene;
        
        %% Debug info below
        if (debug)
            displayDotsPerMeter = displayGet(displayOBJ, 'dots per meter')
            displaySize = displayGet(displayOBJ,'size')
            displayRes = round(displaySize * displayDotsPerMeter)
            imageSize = [size(sensorImageRGB,1) size(sensorImageRGB,1)]
        end
        
        if (debug) && (~isempty(cal))
            xyzImage = sceneGet(scene, 'xyz');
            X1 = squeeze(sensorImageXYZ(:,:,1));
            Y1 = squeeze(sensorImageXYZ(:,:,2));
            Z1 = squeeze(sensorImageXYZ(:,:,3));
            X2 = squeeze(xyzImage(:,:,1));
            Y2 = squeeze(xyzImage(:,:,2));
            Z2 = squeeze(xyzImage(:,:,3));
            [max(X1(:)) max(X2(:)) max(Y1(:)) max(Y2(:)) max(Z1(:)) max(Z2(:))]
        end
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
