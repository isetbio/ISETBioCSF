function run_MosaicsVaryConditionsReviewer
% This is the script used to assess the impact of different mosaic models on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = true;
    makeMosaicsFigure = ~true;
    
    % Mosaic to use
    examinedMosaicModels = {...
        'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection' ...
        'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection' ...
        };
    

    examinedMosaicLegends = {...
        '250,000 cones/mm2' ...
        '288,675 cones/mm2' ...
    };

    idx = [1 2];
    examinedMosaicModels = {examinedMosaicModels{idx}};
     
    % Tun the mosaic-vary condition using the Geisler optics
    opticsName = 'Geisler';
      
    theMosaicTypesAtSpecificSF = {};
    
    % Go !
    for mosaicIndex = 1:numel(examinedMosaicModels)
        mosaicName = examinedMosaicModels{mosaicIndex};
        params = getCSFpaperDefaultParams(mosaicName, computationInstance);
        
        if (strcmp(examinedMosaicLegends{mosaicIndex}, '288,675 cones/mm2'))
            % Denote new mosaics with higher density
            params.mosaicRotationDegs = 360;
        end
        params.opticsModel = opticsName;
       
        params.coneContrastDirection = 'L+M+S';
        params.cyclesPerDegreeExamined = [2 4 8 16 32 50 60];
    
        % Response duration params
        params.frameRate = 10; %(1 frames)
        params.responseStabilizationMilliseconds = 40;
        params.responseExtinctionMilliseconds = 40;

        % Eye movement params
        params.emPathType = 'frozen0';
        params.centeredEMpaths = ~true;

        % Simulation steps to perform
        params.computeMosaic = ~true; 
        params.visualizeMosaic = ~true;

        params.computeResponses = ~true;
        params.computePhotocurrentResponseInstances = ~true;
        params.visualizeMosaic = makeMosaicsFigure;
        params.visualizeResponses = ~true;
        params.visualizeSpatialScheme = ~true;
        params.visualizeOIsequence = ~true;
        params.visualizeOptics = ~true;
        params.visualizeMosaicWithFirstEMpath = ~true;
        params.visualizeSpatialPoolingScheme = ~true;
        params.visualizeStimulusAndOpticalImage = ~true;
        params.visualizeDisplay = ~true;
    
        params.visualizeKernelTransformedSignals = ~true;
        params.findPerformance = ~true;
        params.visualizePerformance = makeSummaryFigure;
        params.deleteResponseInstances = ~true;

        [theMosaicTypes, thePsychometricFunctions{mosaicIndex}, theFigData{mosaicIndex}] = ...
            run_BanksPhotocurrentEyeMovementConditions(params);
        
        if (makeMosaicsFigure)
            for sfIndex = 1:numel(params.cyclesPerDegreeExamined)
                theCurrentMosaic = theMosaicTypes.theMosaics{sfIndex};
                [innerSegmentCoverage(mosaicIndex,sfIndex), geometricCoverage(mosaicIndex,sfIndex)] = ...
                    theCurrentMosaic.retinalCoverage();
            end
            sfIndex = 1;
            theCurrentMosaic = theMosaicTypes.theMosaics{sfIndex};
            theMosaicTypesAtSpecificSF{numel(theMosaicTypesAtSpecificSF) + 1} = theCurrentMosaic;
            theCurrentMosaic.displayInfo();
        end
    end
    
    if (makeSummaryFigure)
        variedParamName = 'MosaicVaryingConeDensity';
        theRatioLims = [0.75 1.45];
        theRatioTicks = 0.8 : 0.1 : 1.4;
        generateFigureForPaper(theFigData, examinedMosaicLegends, variedParamName, opticsName, ...
            'figureType', 'CSF', ...
            'inGraphText', ' E ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end

    if (makeMosaicsFigure)
        hFig = figure(99); clf;
        set(hFig, 'Color', [1 1 1]);
        faceColors = [0.5 0.5 0.5; 1 0.5 0.5];
        edgeColors = [0 0 0; 1 0 0];
        for mosaicIndex = 1:numel(examinedMosaicModels)
            coverage = squeeze(innerSegmentCoverage(mosaicIndex,:));
            normalizedCoverage = coverage / max(coverage);
            subplot(1,2,1)
            hold on;
            plot(params.cyclesPerDegreeExamined, coverage, ...
                'ko-', 'MarkerSize', 14, ...
                'MarkerFaceColor', squeeze(faceColors(mosaicIndex,:)),...
                'MarkerEdgeColor', squeeze(edgeColors(mosaicIndex,:)));
            subplot(1,2,2)
            hold on;
            plot(params.cyclesPerDegreeExamined, normalizedCoverage, ...
                'ko-', 'MarkerSize', 14, ...
                'MarkerFaceColor', squeeze(faceColors(mosaicIndex,:)),...
                'MarkerEdgeColor', squeeze(edgeColors(mosaicIndex,:)));
            
        end % mosaicIdex
        subplot(1,2,1)
        set(gca, 'FontSize', 18, 'XLim', [0.5 90], 'YLim', [0 1]);
        hL = legend(examinedMosaicLegends);
        grid on;
        box on;
        xlabel('\it spatial frequency (c/deg)');
        ylabel('\it coverage factor');
        
        subplot(1,2,2)
        set(gca, 'FontSize', 18, 'XLim', [0.5 90], 'YLim', [0.7 1]);
        hL = legend(examinedMosaicLegends);
        grid on;
        box on;
        xlabel('\it spatial frequency (c/deg)');
        ylabel('\it normalized coverage factor');
        drawnow
        
        generateMosaicsFigure(theMosaicTypesAtSpecificSF, examinedMosaicLegends,  ...
            'inGraphTexts', {' A ', ' B ', ' C ', ' D '}, ...
            'visualizedFOV', 0.35);
        
    end
end

function generateMosaicsFigure(theMosaics, examinedMosaicLegends,  varargin)
    p = inputParser;
    p.addParameter('inGraphTexts', '', @iscell);
    p.addParameter('visualizedFOV', 0.5, @isnumeric);
    p.parse(varargin{:});
    
    inGraphTexts = p.Results.inGraphTexts;
    visualizedFOV = p.Results.visualizedFOV;
    
    % Initialize figure
    hFig = figure(1); clf;
    formatFigureForPaper(hFig, 'figureType', 'MOSAICS');
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 2, ...
       'heightMargin',  0.01, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.01, ...
       'rightMargin',    0.0001, ...
       'bottomMargin',   0.06, ...
       'topMargin',      0.0001);
    
    for mosaicIndex = 1:numel(theMosaics)
        cm = theMosaics{mosaicIndex};
        posRangeY = 0.7*visualizedFOV*cm.micronsPerDegree * 1e-6 * [-1 1];
        posRangeX = posRangeY;
        mosaicName = examinedMosaicLegends{mosaicIndex};
        
        switch (mosaicIndex)
            case 1
                subplotRow = 1; subplotCol = 1;
            case 2
                subplotRow = 1; subplotCol = 2;
            case 3
                subplotRow = 2; subplotCol = 1;
            case 4
                subplotRow = 2; subplotCol = 2;
        end
        
        ax = subplot('Position', subplotPosVectors(subplotRow,subplotCol).v);
        cm.visualizeGrid(...
            'axesHandle', ax, ...
            'apertureShape', 'disks', ...
            'visualizedConeAperture', 'lightCollectingArea', ...
            'labelConeTypes', true, ...
            'overlayHexMesh', false, ...
            'backgroundcolor', [0 0 0]);
        
        set(ax, 'XLim', posRangeX, 'YLim', posRangeY);
        
        inGraphTextPos = [-0.25 0.21]*cm.micronsPerDegree * 1e-6;
        t = text(ax, inGraphTextPos(1), inGraphTextPos(2), inGraphTexts{mosaicIndex});
        
        formatFigureForPaper(hFig, ...
            'figureType', 'MOSAICS', ...
            'theAxes', ax, ...
            'theText', t, ...
            'theTextFontSize', [], ...
            'theFigureTitle', '');
        
        tickDegs = -0.5:0.1:0.5;
        ticks = tickDegs*cm.micronsPerDegree * 1e-6;
        tickLabels = sprintf('%2.2f\n', tickDegs);
        set(ax, 'XTick', ticks, 'XTickLabel', tickLabels, 'YTick', ticks, 'YTickLabel', tickLabels);
        
        switch (mosaicIndex)
            case 1
                set(ax, 'XTickLabel', {});
                set(ax, 'YTickLabel', {});
            case 2
                set(ax, 'XTickLabel', {});
                set(ax, 'YTickLabel', {});
            case 3
                xlabel(ax, 'retinal position (deg)', 'FontWeight', 'bold');
                set(ax, 'YTickLabel', {});
            case 4
                xlabel(ax, 'retinal position (deg)', 'FontWeight', 'bold');
                set(ax, 'YTickLabel', {});
        end
    end
    
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    figureName = fullfile(exportsDir, sprintf('MosaicsEmployed.pdf'));
    NicePlot.exportFigToPDF(figureName, hFig, 300);
end

