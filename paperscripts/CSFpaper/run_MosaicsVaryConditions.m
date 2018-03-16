function run_MosaicsVaryConditions
% This is the script used to assess the impact of different mosaic models on the CSF
%  
    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    % Whether to make a summary figure with CSF from all examined conditions
    makeSummaryFigure = ~true;
    makeMosaicsFigure = true;
    
    % Mosaic to use
    examinedMosaicModels = {...
        'originalBanks' ...
        'ISETbioHexEccBasedNoScones' ...
        'ISETbioHexEccBasedLMSrealistic' ...
        };
 
    examinedMosaicLegends = {...
        'constant LM density (Banks ''87)' ...
        'eccentricity-based LM density' ...
        'eccentricity-based LMS density' ...
    };

    % Tun the mosaic-vary condition using the Geisler optics
    opticsName = 'Geisler';
      
    % Go !
    for mosaicIndex = 1:numel(examinedMosaicModels)
        mosaicName = examinedMosaicModels{mosaicIndex};
        params = getCSFpaperDefaultParams(mosaicName, computationInstance);
        
        % Special case: the Geisler optics/original Banks mosaic was run up
        % to 50 c/deg
        if (strcmp(mosaicName, 'originalBanks')) && (strcmp(opticsName,'Geisler'))
            params.cyclesPerDegreeExamined = params.cyclesPerDegreeExamined(1:end-1);
        end
        
        params.opticsModel = opticsName;
       
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
    
        params.visualizeKernelTransformedSignals = ~true;
        params.findPerformance = ~true;
        params.visualizePerformance = makeSummaryFigure;
        params.deleteResponseInstances = ~true;

        [theMosaics{mosaicIndex},thePsychometricFunctions{mosaicIndex}, theFigData{mosaicIndex}] = ...
            run_BanksPhotocurrentEyeMovementConditions(params);
    end
    
    if (makeSummaryFigure)
        variedParamName = 'Mosaic';
        theRatioLims = [0.3 1.2];
        theRatioTicks = [0.3 0.5 0.7 1.0 2.0];
        generateFigureForPaper(theFigData, examinedMosaicLegends, variedParamName, opticsName, ...
            'figureType', 'CSF', ...
            'inGraphText', ' A ', ...
            'plotFirstConditionInGray', true, ...
            'plotRatiosOfOtherConditionsToFirst', true, ...
            'theRatioLims', theRatioLims, ...
            'theRatioTicks', theRatioTicks ...
            );
    end
    
    if (makeMosaicsFigure)
        generateMosaicsFigure(theMosaics, examinedMosaicLegends,  ...
            'inGraphTexts', {'B', 'C', 'D'}, ...
            'visualizedFOV', 0.3);
    end
end

function generateMosaicsFigure(theMosaics, examinedMosaicLegends,  varargin)
    p = inputParser;
    p.addParameter('inGraphTexts', '', @iscell);
    p.addParameter('visualizedFOV', 0.4, @isnumeric);
    p.parse(varargin{:});
    
    inGraphTexts = p.Results.inGraphTexts;
    visualizedFOV = p.Results.visualizedFOV;
    
    % Initialize figure
    hFig = figure(1); clf;
    formatFigureForPaper(hFig, 'figureType', 'MOSAICS');
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', numel(theMosaics), ...
       'colsNum', 1, ...
       'heightMargin',  0.05, ...
       'widthMargin',    0.001, ...
       'leftMargin',     0.01, ...
       'rightMargin',    0.001, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.01);
    
    for mosaicIndex = 1:numel(theMosaics)
        cm = theMosaics{mosaicIndex};
        
        posRangeY = 0.5*visualizedFOV*cm.micronsPerDegree * 1e-6 * [-1 1];
        posRangeX = 2*posRangeY;
        mosaicName = examinedMosaicLegends{mosaicIndex};
        ax = subplot('Position', subplotPosVectors(mosaicIndex,1).v);
        cm.visualizeGrid(...
            'axes handle', ax, ...
            'aperture shape', 'disks', ...
            'visualizedConeAperture', 'geometricArea', ...
            'label cone types', true, ...
            'overlayHexMesh', false, ...
            'backgroundcolor', [1 1 1]);
        
        formatFigureForPaper(hFig, ...
            'figureType', 'MOSAICS', ...
            'inGraphText', inGraphTexts{mosaicIndex}, ...
            'theAxes', ax, ...
            'theFigureTitle', mosaicName);
        
        set(ax, 'XLim', posRangeX, 'YLim', posRangeY);
        tickDegs = -0.5:0.1:0.5
        ticks = tickDegs*cm.micronsPerDegree * 1e-6
        tickLabels = sprintf('%2.2f\n', ticks)
        set(ax, 'XTick', ticks, 'XTickLabel', tickLabels, 'YTick', ticks, 'YTickLabel', tickLabels);
        if (mosaicIndex == numel(theMosaics))
            xlabel(ax, 'deg');
        else
            set(ax, 'XTickLabel', {});
        end
        ylabel(ax, 'deg');
        
        title(ax,mosaicName);
    end
    
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    figureName = fullfile(exportsDir, sprintf('MosaicsEmployed.pdf'));
    NicePlot.exportFigToPDF(figureName, hFig, 300);
end

