function analyzePowerSpectrumChangeDueToDrift

    figNo = 0;
    
    
    load('1overF.mat', 'oneOverFnoiseFEM', 'oneOverFnoiseStabilized');
    figNo = figNo + 1;
    analyzeSpectra(oneOverFnoiseFEM, oneOverFnoiseStabilized, figNo);
%     
%     load('1overF_smallerWindows.mat', 'oneOverFnoiseFEM', 'oneOverFnoiseStabilized');
%     figNo = figNo + 1;
%     analyzeSpectra(oneOverFnoiseFEM, oneOverFnoiseStabilized, figNo);
%     
%     
%     load('lowFrequency.mat', 'lowFreqFEM', 'lowFreqStabilized');
%     figNo = figNo + 1;
%     analyzeSpectra(lowFreqFEM, lowFreqStabilized, figNo);
%     
%     load('lowFrequency_smallerWindows.mat', 'lowFreqFEM', 'lowFreqStabilized');
%     figNo = figNo + 1;
%     analyzeSpectra(lowFreqFEM, lowFreqStabilized, figNo);
        
    load('highFrequency.mat', 'highFreqFEM', 'highFreqStabilized');
    figNo = figNo + 1;
    analyzeSpectra(highFreqFEM, highFreqStabilized, figNo);
    
%     load('highFrequency_smallerWindows.mat', 'highFreqFEM', 'highFreqStabilized');
%     figNo = figNo + 1;
%     analyzeSpectra(highFreqFEM, highFreqStabilized, figNo);
    
    
end

