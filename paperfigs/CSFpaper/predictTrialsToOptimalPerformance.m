function predictTrialsToOptimalPerformance
    
    
    figure(1); clf;
    
    load('/Users/nicolas/Desktop/theData8.mat', 'theData')
    subplot(1,3,1); 
    mlptData = theData{1};
    svmPCA = theData{2};
    svmPool = theData{3};
    trialsLog = log10(mlptData.trialsNum);
    mlpt.csLog = log10(mlptData.contrastSensitivity);
    svmPCA.csLog = log10(svmPCA.contrastSensitivity);
    svmPool.csLog = log10(svmPool.contrastSensitivity);
    
    coeffs = polyfit(trialsLog, svmPCA.csLog, 1)
    % Get fitted values
    svmPCA.csLogExtrapX = linspace(2, 30, 500);
    svmPCA.csLogExtrapY = polyval(coeffs, svmPCA.csLogExtrapX);
    % Get fitted values
    coeffs = polyfit(trialsLog, svmPool.csLog, 1);
    svmPool.csLogExtrapX = svmPCA.csLogExtrapX
    svmPool.csLogExtrapY = polyval(coeffs, svmPool.csLogExtrapX);
    
    plot([2 30], mlpt.csLog(1)*[1 1], 'k--'); hold on;
    plot(trialsLog, svmPCA.csLog, 'ks-'); hold on;
    plot(trialsLog, svmPool.csLog, 'rs-'); hold on;
    plot(svmPCA.csLogExtrapX, svmPCA.csLogExtrapY, 'b.-');
    plot(svmPool.csLogExtrapX, svmPool.csLogExtrapY, 'r.-');
    
    
    set(gca, 'XLim', [2 20], 'YLim', [0.1 3]);
    
    load('/Users/nicolas/Desktop/theData16.mat', 'theData')
    subplot(1,3,2);
    mlptData = theData{1};
    svmPCA = theData{2};
    svmPool = theData{3};
    trialsLog = log10(mlptData.trialsNum);
    mlpt.csLog = log10(mlptData.contrastSensitivity);
    svmPCA.csLog = log10(svmPCA.contrastSensitivity);
    svmPool.csLog = log10(svmPool.contrastSensitivity);
    coeffs = polyfit(trialsLog, svmPCA.csLog, 1)
    % Get fitted values
    svmPCA.csLogExtrapX = linspace(2, 30, 500);
    svmPCA.csLogExtrapY = polyval(coeffs, svmPCA.csLogExtrapX);
    coeffs = polyfit(trialsLog, svmPool.csLog, 1);
    % Get fitted values
    svmPool.csLogExtrapX = svmPCA.csLogExtrapX
    svmPool.csLogExtrapY = polyval(coeffs, svmPool.csLogExtrapX);
    
    
    plot([2 30], mlpt.csLog(1)*[1 1], 'k--'); hold on;
    plot(trialsLog, svmPCA.csLog, 'ks-'); hold on;
    plot(trialsLog, svmPool.csLog, 'rs-'); hold on;
    plot(svmPCA.csLogExtrapX, svmPCA.csLogExtrapY, 'b.-');
    plot(svmPool.csLogExtrapX, svmPool.csLogExtrapY, 'r.-');
    set(gca, 'XLim', [2 20], 'YLim', [0.1 3]);
    
    load('/Users/nicolas/Desktop/theData32.mat', 'theData')
    subplot(1,3,3);
    mlptData = theData{1};
    svmPCA = theData{2};
    svmPool = theData{3};
    trialsLog = log10(mlptData.trialsNum);
    mlpt.csLog = log10(mlptData.contrastSensitivity);
    svmPCA.csLog = log10(svmPCA.contrastSensitivity);
    svmPool.csLog = log10(svmPool.contrastSensitivity);
    coeffs = polyfit(trialsLog, svmPCA.csLog, 1)
    % Get fitted values
    svmPCA.csLogExtrapX = linspace(2,30, 500);
    svmPCA.csLogExtrapY = polyval(coeffs, svmPCA.csLogExtrapX);
    coeffs = polyfit(trialsLog, svmPool.csLog, 1);
    % Get fitted values
    svmPool.csLogExtrapX = svmPCA.csLogExtrapX
    svmPool.csLogExtrapY = polyval(coeffs, svmPool.csLogExtrapX);
    
    
    plot([2 30], mlpt.csLog(1)*[1 1], 'k--'); hold on;
    plot(trialsLog, svmPCA.csLog, 'ks-'); hold on;
    plot(trialsLog, svmPool.csLog, 'rs-'); hold on;
    plot(svmPCA.csLogExtrapX, svmPCA.csLogExtrapY, 'b.-');
    plot(svmPool.csLogExtrapX, svmPool.csLogExtrapY, 'r.-');
    set(gca, 'XLim', [2 20], 'YLim', [0.1 3]);
end
