function predictTrialsToOptimalPerformance
    
    load('/Users/nicolas/Desktop/theData8.mat', 'theData')
    [trialsLog8, mlpt8, svmPCA8, svmPool8] = analyzeData(theData);
    
    load('/Users/nicolas/Desktop/theData16.mat', 'theData')
    [trialsLog16, mlpt16, svmPCA16, svmPool16] = analyzeData(theData);
    
    load('/Users/nicolas/Desktop/theData32.mat', 'theData')
    [trialsLog32, mlpt32, svmPCA32, svmPool32] = analyzeData(theData);
    
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1000 1000]);
    subplot(2,2,1); 
    plotData(trialsLog8, mlpt8, svmPCA8, svmPool8, '8 c/deg');
    
    subplot(2,2,2); 
    plotData(trialsLog16, mlpt16, svmPCA16, svmPool16, '16 c/deg');
    
    subplot(2,2,3); 
    plotData(trialsLog32, mlpt32, svmPCA32, svmPool32, '32 c/deg');
    
end

function [trialsLog, mlpt, svmPCA, svmPool] = analyzeData(theData)
    mlptData = theData{1};
    svmPCA = theData{2};
    svmPool = theData{3};
    trialsLog = log10(mlptData.trialsNum);
    mlpt.csLog = log10(mlptData.contrastSensitivity);
    svmPCA.csLog = log10(svmPCA.contrastSensitivity);
    svmPool.csLog = log10(svmPool.contrastSensitivity);
    
    coeffs = polyfit(trialsLog, svmPCA.csLog, 1);
    svmPCA.csLogExtrapX = linspace(2, 30, 500);
    svmPCA.csLogExtrapY = polyval(coeffs, svmPCA.csLogExtrapX);
    
    coeffs = polyfit(trialsLog, svmPool.csLog, 1);
    svmPool.csLogExtrapX = svmPCA.csLogExtrapX;
    svmPool.csLogExtrapY = polyval(coeffs, svmPool.csLogExtrapX);
    
end

function plotData(trialsLog, mlpt, svmPCA, svmPool, titleText)
    
    plot([2 30], mlpt.csLog(1)*[1 1], 'k--', 'MarkerSize', 12); hold on;
    plot(trialsLog, svmPCA.csLog, 'ks-', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 0.50]); hold on;
    plot(trialsLog, svmPool.csLog, 'rs-', 'MarkerSize', 12, 'MarkerFaceColor', [0.9 0.5 0.50]); hold on;
    plot(svmPCA.csLogExtrapX, svmPCA.csLogExtrapY, 'k-', 'LineWidth', 1.5);
    plot(svmPool.csLogExtrapX, svmPool.csLogExtrapY, 'r-', 'LineWidth', 1.5);
    hL = legend('MLPT', 'SVM-PCA', 'SVM-Pool', 'SVM-PCA fit', 'SVM-Pool fit');
    set(hL, 'Location', 'SouthEast', 'FontSize', 14)
    set(gca, 'XLim', [2 8], 'YLim', [0.1 3], 'FontSize', 16);
    grid on; box on
    xlabel('log10 trials', 'FontWeight', 'bold');
    ylabel('log10 CSF', 'FontWeight', 'bold');
    axis 'square'
    title(titleText);
end
