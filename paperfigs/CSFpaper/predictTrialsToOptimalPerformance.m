function predictTrialsToOptimalPerformance
    
    load('/Users/nicolas/Desktop/theData8.mat', 'theData')
    [trialsLog8, mlpt8, svmPCA8, svmPool8] = analyzeData(theData);
    
    load('/Users/nicolas/Desktop/theData16.mat', 'theData')
    [trialsLog16, mlpt16, svmPCA16, svmPool16] = analyzeData(theData);
    
    load('/Users/nicolas/Desktop/theData32.mat', 'theData')
    [trialsLog32, mlpt32, svmPCA32, svmPool32] = analyzeData(theData);
    
    figure(1); clf;
    subplot(1,3,1); 
    plotData(trialsLog8, mlpt8, svmPCA8, svmPool8);
    
    subplot(1,3,2); 
    plotData(trialsLog16, mlpt16, svmPCA16, svmPool16);
    
    subplot(1,3,3); 
    plotData(trialsLog32, mlpt32, svmPCA32, svmPool32);
    
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

function plotData(trialsLog, mlpt, svmPCA, svmPool)
    
    plot([2 30], mlpt.csLog(1)*[1 1], 'k--'); hold on;
    plot(trialsLog, svmPCA.csLog, 'ks-'); hold on;
    plot(trialsLog, svmPool.csLog, 'rs-'); hold on;
    plot(svmPCA.csLogExtrapX, svmPCA.csLogExtrapY, 'b.-');
    plot(svmPool.csLogExtrapX, svmPool.csLogExtrapY, 'r.-');
    set(gca, 'XLim', [2 20], 'YLim', [0.1 3]);
end
