function importColorMaterialISETbioData
% Loads the isomerization responses to the set of the 13 images employed
% in the the color-material paper and displays them.
%
% History
% 2/8/18 NPC Wrote it.
%

    close all;
    
    % Choose whether to render images of the data or not
    renderIsomerizationMaps = ~true;
    renderDemosaiceIsomerizationMaps = true;
             
    % Load the data
    horizontalFOV = 15;
    allScenesData = loadAllSceneData(horizontalFOV);
    
    % Retrieve the scene names
    scenesList = keys(allScenesData);
    
    % Loop through all scenes
    for sceneIndex = 1:numel(scenesList)
         theSceneName = scenesList{sceneIndex};
         
         % Retrieve the data for this scene
         theSceneData = allScenesData(theSceneName);
         
         % Retrieve the cone isomerization rates and cone locations
         isomerizationRates = theSceneData.theConeIsomerizations;
         xConeLocsDegs = theSceneData.xConeLocsDegs;
         yConeLocsDegs = theSceneData.yConeLocsDegs;
         
         % Retrieve the demosaiced isomerization rates and their spatial support
         demosaicedIsomerizationsMaps = theSceneData.theDemosaicedIsomerizationsMaps;
         demosaicedIsomerizationsMapsSupport = theSceneData.theDemosaicedIsomerizationsMapSupport;
         
         % Retrieve the max isomerization rate for this scene
         maxIsomerizationRate = theSceneData.maxIsomerizationRate;
         
         % Retrieve the cone names
         coneTypeNames = keys(isomerizationRates);
         
         % Figure handles
         hFig1 = []; hFig2 = [];
         
         % Loop throught the cone types, retrieve the data and plot it
         for coneTypeIndex = 1:numel(coneTypeNames)
             % Retrieve the cone type name
             coneTypeName = coneTypeNames{coneTypeIndex};
             
             % Retrieve the isomerization rates for this cone type
             isomerizationRatesForThisConeType = isomerizationRates(coneTypeName);
             xConeLocsDegsForThisConeType = xConeLocsDegs(coneTypeName);
             yConeLocsDegsForThisConeType = yConeLocsDegs(coneTypeName);
             
             if (sceneIndex == 1)
                 fprintf('%s: %d\n', coneTypeName, numel(isomerizationRatesForThisConeType));
             end
             
             % Retrieve the demosaiced isomerization map for this cone type
             demosaicedIsomerizationsMapForThisConeType = demosaicedIsomerizationsMaps(coneTypeName);
             demosaicedIsomerizationsMapSupportForThisConeType = demosaicedIsomerizationsMapsSupport(coneTypeName);
             
             if (renderIsomerizationMaps)
                 % Plot the isomerization maps
                 if (sceneIndex == 1) && (coneTypeIndex == 1)
                    coneApertureOutline = generateConeApertureOutline(horizontalFOV);
                 end
                 hFig1 = plotMosaicActivationMaps(isomerizationRatesForThisConeType, ...
                     xConeLocsDegsForThisConeType, yConeLocsDegsForThisConeType, maxIsomerizationRate, ...
                     coneTypeName, coneTypeIndex, coneApertureOutline, hFig1, theSceneName);
                 NicePlot.exportFigToPNG(sprintf('../%s_Isomerizations.png',theSceneName), hFig1, 300);
                 %close(hFig1);
             end
             
             if (renderDemosaiceIsomerizationMaps)
                 % Plot the demosaiced isomerization maps
                 hFig2 = plotDemosaicedActivationMaps(demosaicedIsomerizationsMapForThisConeType, ...
                     demosaicedIsomerizationsMapSupportForThisConeType, maxIsomerizationRate, ...
                     coneTypeName, coneTypeIndex, hFig2, theSceneName);
                 NicePlot.exportFigToPNG(sprintf('../%s_DemosaicedIsomerizations.png',theSceneName), hFig2, 300);
                 %close(hFig2);
             end
         end % coneTypeIndex
    end % sceneIndex
end

% Method to plot the isomeration responses
function hFig = plotMosaicActivationMaps(isomerizationRates, xConeLocsDegs, yConeLocsDegs, maxIsomerizationRate, coneTypeName, coneTypeIndex, coneApertureOutline, hFig, theSceneName)
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 4, ...
       'heightMargin',   0.06, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.006, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);
    
    if isempty(hFig)
        hFig = figure(); clf;
        set(hFig, 'Position', [10 10 1670 920], 'Name', theSceneName);
        colormap(bone(1024));
    end
               
    figure(hFig);
    for normalizationModeIndex = 1:2
        subplot('Position', subplotPosVectors(normalizationModeIndex,coneTypeIndex).v);
        if (normalizationModeIndex == 1)
            renderPatchArray(gca, coneApertureOutline, xConeLocsDegs, yConeLocsDegs, isomerizationRates/maxIsomerizationRate, 'none', 0.5);
        else
            renderPatchArray(gca, coneApertureOutline, xConeLocsDegs, yConeLocsDegs, isomerizationRates/max(isomerizationRates(:)), 'none', 0.5);
        end
        axis 'image'; axis 'xy';
        mD = max([max(abs(xConeLocsDegs)) max(abs(yConeLocsDegs))]);
        if (mD > 3)
            ticks = -15:1:15;
        else
            ticks = -15:0.25:15;
        end
        range = 1.05 * mD * [-1 1];
        set(gca, 'CLim', [0 1], 'XTick', ticks, 'YTick', ticks, ...
                 'XLim', range, 'YLim', range, 'FontSize', 12);
        box on; grid on;
        xlabel('space (degs');
        title(sprintf('%s (max: %2.0f R*/sec)', coneTypeName, max(isomerizationRates(:))));
        if (coneTypeIndex>1)
            set(gca, 'YTickLabel',{});
        end
    end % normalizationModeIndex
    
    % Nested function for efficiently rendering large number of patches
    function renderPatchArray(axesHandle, pixelOutline, xCoords, yCoords, faceColorsNormalizedValues,  edgeColor, lineWidth)
        verticesPerCone = numel(pixelOutline.x);
        verticesList = zeros(verticesPerCone * numel(xCoords), 2);
        facesList = []; colors = [];
        for coneIndex = 1:numel(xCoords)
            idx = (coneIndex-1)*verticesPerCone + (1:verticesPerCone);
            verticesList(idx,1) = pixelOutline.x(:) + xCoords(coneIndex);
            verticesList(idx,2) = pixelOutline.y(:) + yCoords(coneIndex);
            facesList = cat(1, facesList, idx);
            colors = cat(1, colors, repmat(faceColorsNormalizedValues(coneIndex), [verticesPerCone 1]));
        end

        S.Vertices = verticesList;
        S.Faces = facesList;
        S.FaceVertexCData = colors;
        S.FaceColor = 'flat';
        S.EdgeColor = edgeColor;
        S.LineWidth = lineWidth;
        patch(S, 'Parent', axesHandle);
    end
end

% Method to plot the demosaiced isomerization map
function hFig = plotDemosaicedActivationMaps(demosaicedMap, demosaicedMapSupport, maxMapActivation, coneTypeName, coneTypeIndex, hFig, theSceneName)
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 4, ...
       'heightMargin',   0.06, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.006, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);
    
    if isempty(hFig)
        hFig = figure(); clf;
        set(hFig, 'Position', [10 10 1670 920], 'Name', theSceneName);
        colormap(bone(1024));
    end
         
    figure(hFig);
    for normalizationModeIndex = 1:2
        subplot('Position', subplotPosVectors(normalizationModeIndex,coneTypeIndex).v);
        if (normalizationModeIndex == 1)
            imagesc(demosaicedMapSupport, demosaicedMapSupport, demosaicedMap/maxMapActivation);
        else
            imagesc(demosaicedMapSupport, demosaicedMapSupport, demosaicedMap/max(demosaicedMap(:)));
        end
        axis 'image'; axis 'xy';
        if (max(demosaicedMapSupport) > 3)
            ticks = -15:1:15;
        else
            ticks = -15:0.25:15;
        end
        range = 1.05*[demosaicedMapSupport(1) demosaicedMapSupport(end)];
        set(gca, 'CLim', [0 1], 'XTick', ticks, 'YTick', ticks, ...
                 'XLim', range, 'YLim', range, 'FontSize', 12);
        box on; grid on;
        xlabel('space (degs');
        title(sprintf('%s (max: %2.0f R*/sec)', coneTypeName, max(demosaicedMap(:))));
        if (coneTypeIndex>1)
            set(gca, 'YTickLabel',{});
        end
    end                     
end

% Method to load all the data
function allScenesData = loadAllSceneData(horizontalFOV)
    responsesMatFileName = sprintf('../responses_%2.2fdegFOV.mat', horizontalFOV);
    fprintf('\n\nLoading responses from %s ...', responsesMatFileName);
    load(responsesMatFileName, 'allScenesData');
    fprintf('Finished importing responses.\n');
end

function apertureOutline = generateConeApertureOutline(horizontalFOV)
    coneMosaicMatFileName = sprintf('../coneMosaic_%2.2fdegFOV.mat', horizontalFOV);
    load(coneMosaicMatFileName, 'theConeMosaic');
    coneApertureDiameterMicrons = diameterForCircularApertureFromWidthForSquareAperture(theConeMosaic.pigment.pdWidth*1e6);
    
    radiusDeg = 0.5*coneApertureDiameterMicrons/theConeMosaic.micronsPerDegree;
                
    iTheta = (0:60:360)/180*pi;
    apertureOutline.x = radiusDeg * cos(iTheta);
    apertureOutline.y = radiusDeg * sin(iTheta);
end
