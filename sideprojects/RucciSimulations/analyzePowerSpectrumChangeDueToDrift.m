function analyzePowerSpectrumChangeDueToDrift

    load('1overF.mat', 'oneOverFnoiseStimulus', 'oneOverFnoiseFEM', 'oneOverFnoiseStabilized');
    figNo = 1;
    analyzeSpectra(oneOverFnoiseFEM, oneOverFnoiseStabilized, figNo);
    
end

