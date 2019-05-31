function analyzePowerSpectrumChangeDueToDrift

    load('1overF.mat', 'oneOverFnoiseFEM', 'oneOverFnoiseStabilized');
    figNo = 1;
    analyzeSpectra(oneOverFnoiseFEM, oneOverFnoiseStabilized, figNo);
    
    load('1overF_smallerWindows.mat', 'oneOverFnoiseFEM', 'oneOverFnoiseStabilized');
    figNo = 2;
    analyzeSpectra(oneOverFnoiseFEM, oneOverFnoiseStabilized, figNo);
    
    
%     load('lowFrequency.mat', 'lowFrequencyStimulus', 'lowFreqFEM', 'lowFreqStabilized');
%     figNo = figNo + 1;
%     analyzeSpectra(lowFreqFEM, lowFreqStabilized, figNo);
%     
%     load('highFrequency.mat', 'highFrequencyStimulus', 'highFreqFEM', 'highFreqStabilized');
%     figNo = figNo + 1;
%     analyzeSpectra(highFreqFEM, highFreqStabilized, figNo);
    
    
end

