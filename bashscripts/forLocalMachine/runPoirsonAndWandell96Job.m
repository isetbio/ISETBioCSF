function runPoirsonAndWandell96Job()

    
    %classifierTypeList = {'svm'};
    pcaComponentsNum = 60;
    
    
    % Actions to perform
    computeMosaic = false;
    computeResponses   = true;
    findPerformances   = true;
    
    % Visualization options
    visualizeSpatialScheme = ~true;
    visualizeResponses = ~true;
    visualizePerformances = true;
    visualizeMosaic = ~true;
    
    nTrials = 1024*4;
    performanceClassifierTrainingSamples = 1024*4;
    
    
    % Redo findPerformance (first for energy, then for half-wave-rect), then use below comparison plots
    for configID = [2 5]
    
        classifierSignalList = {'isomerizations', 'photocurrents'};
        classifierTypeList = {'svmV1FilterBank'};
        
        c_PoirsonAndWandell96RunSession(configID, nTrials, performanceClassifierTrainingSamples, ...
            computeMosaic, computeResponses, findPerformances, visualizeResponses, visualizePerformances, ...
            visualizeMosaic, visualizeSpatialScheme, classifierSignalList, classifierTypeList, pcaComponentsNum);

    end
    return;
    
    
    % Compare emRandom vs emFrozen0
    configID1 = 4; configID2 = 1;
    %configID1 = 5; configID2 = 2;
    %configID1 = 6; configID2 = 3;
    
    classifierSignalList = {'isomerizations'};
    classifierTypeList = {'mlpt'};
    d1 = c_PoirsonAndWandell96RunSession(configID1, nTrials, performanceClassifierTrainingSamples, ...
        computeMosaic, computeResponses, findPerformances, visualizeResponses, visualizePerformances, ...
        visualizeMosaic, visualizeSpatialScheme, classifierSignalList, classifierTypeList, pcaComponentsNum);
    
    classifierSignalList = {'isomerizations', 'photocurrents'};
    classifierTypeList = {'svmV1FilterBank'};
    
	d2 = c_PoirsonAndWandell96RunSession(configID1, nTrials, performanceClassifierTrainingSamples, ...
        computeMosaic, computeResponses, findPerformances, visualizeResponses, visualizePerformances, ...
        visualizeMosaic, visualizeSpatialScheme, classifierSignalList, classifierTypeList, pcaComponentsNum);
    
    
    
    classifierSignalList = {'isomerizations'};
    classifierTypeList = {'mlpt'};
    d3 = c_PoirsonAndWandell96RunSession(configID2, nTrials, performanceClassifierTrainingSamples, ...
        computeMosaic, computeResponses, findPerformances, visualizeResponses, visualizePerformances, ...
        visualizeMosaic, visualizeSpatialScheme, classifierSignalList, classifierTypeList, pcaComponentsNum);
    
    classifierSignalList = {'isomerizations', 'photocurrents'};
    classifierTypeList = {'svmV1FilterBank'};
    
	d4 = c_PoirsonAndWandell96RunSession(configID2, nTrials, performanceClassifierTrainingSamples, ...
        computeMosaic, computeResponses, findPerformances, visualizeResponses, visualizePerformances, ...
        visualizeMosaic, visualizeSpatialScheme, classifierSignalList, classifierTypeList, pcaComponentsNum);
    
    thresholds = [...
        d1{1}.detectionThresholdContrast d3{1}.detectionThresholdContrast;
        d2{1}.detectionThresholdContrast d4{1}.detectionThresholdContrast;
        d2{2}.detectionThresholdContrast d4{2}.detectionThresholdContrast ...
        ];
    
    
    conditions = {...
        sprintf('%s-%s', d1{1}.performanceSignalName, d1{1}.classifierTypeName) ...
        sprintf('%s-%s', d2{1}.performanceSignalName, d2{1}.classifierTypeName) ...
        sprintf('%s-%s', d2{2}.performanceSignalName, d2{2}.classifierTypeName)  };


    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1, ...
           'colsNum', 1, ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.07, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.02);
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 60 860 1290], 'Color', [1 1 1]);
    
    subplot('Position', subplotPosVectors(1,1).v);
    bPlot = bar([1 2 3], thresholds, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [.3 .3 .3]);
    bPlot(1).FaceColor = [0.8 0.8 0.8];
    bPlot(2).FaceColor = [1.0 0.0 0.0];
    set(gca, 'XTickLabels', conditions, 'FontSize', 16, 'YLim', [0 0.072]);
    hL = legend({'em: frozen0', 'em: random'});
    set(hL, 'FontSize', 14, 'Location', 'NorthWest');
    grid on;
    ylabel('detection threshold contrast', 'FontWeight', 'bold');
    xlabel('signal - classifier', 'FontWeight', 'bold');
    title(sprintf('%2.1f c/deg, %2.1f cd/m2', d1{1}.stimParams.spatialFrequency, d1{1}.stimParams.meanLuminance));
    
    drawnow;
    NicePlot.exportFigToPDF('cond1.pdf', hFig, 300);
end


