function simulatePowerSpectrumChangeDueToDrift

    % Experiment params
    nTrials = 16;
    emDurationSeconds = 1.0;
    
    % Generate fixational eye movmeents
    [emPosArcMin, timeAxis] = generateFEM(nTrials, emDurationSeconds);
    
    % Stimulus spatial params
    coneApertureMicrons = 2; micronsPerDegree = 300;
    stimulusPixelSizeArcMin = 1*coneApertureMicrons / micronsPerDegree * 60;
    stimulusWidthArcMin = 200;
    
    % Grating params
    gratingParams.oriDegs = 90;
    gratingParams.sigmaArcMin = 30;
    gratingParams.contrast = 1;
    gratingParams.sfCPD = 1;
    
    % Noise params
    noiseParams.spectrumShape = '';
    noiseParams.cornerFrequencyCPD = [];
    noiseParams.steepness = 100;

    
    figNo = 1;
    noiseNorm = nan;
    
    
    do_1overF_simulation = ~true;
    do_lowFrequency_simulation = ~true;
    do_highFrequency_simulation = true;
    
    if (do_1overF_simulation)
        % Noise params for noise only stimulus
        gratingParams.contrast = 0.0;
        noiseParams.spectrumShape = '1overF';
    
        [oneOverFnoiseStimulus, noiseNorm] = generateStimulusImage(stimulusWidthArcMin, stimulusPixelSizeArcMin,  gratingParams, noiseParams, noiseNorm, nTrials, figNo);
        oneOverFnoiseFEM = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, oneOverFnoiseStimulus);
        oneOverFnoiseStabilized = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin*0, oneOverFnoiseStimulus);
        save('1overF.mat', 'oneOverFnoiseStimulus', 'oneOverFnoiseFEM', 'oneOverFnoiseStabilized', '-v7.3');
        figNo = figNo + 1;
        analyzeSpectra(oneOverFnoiseFEM, oneOverFnoiseStabilized, figNo);
    end
    
    if (do_lowFrequency_simulation)
        % Params for low frequency stimulus
        gratingParams.sfCPD = 4;
        gratingParams.contrast = 1;
        noiseParams.spectrumShape = 'highPassCornerFreq';
        noiseParams.cornerFrequencyCPD = 10;
    
        [lowFrequencyStimulus, noiseNorm] = generateStimulusImage(stimulusWidthArcMin, stimulusPixelSizeArcMin,  gratingParams, noiseParams, noiseNorm, nTrials, figNo);
        lowFreqFEM = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, lowFrequencyStimulus);
        lowFreqStabilized = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin*0, lowFrequencyStimulus);
        save('lowFrequency.mat', 'lowFrequencyStimulus', 'lowFreqFEM', 'lowFreqStabilized', '-v7.3');
        
        figNo = figNo + 1;
        analyzeSpectra(lowFreqFEM, lowFreqStabilized, figNo);
    end
    
    
    if (do_highFrequency_simulation)
        % Noise params for highfrequency stimulus
        gratingParams.sfCPD = 11;
        gratingParams.contrast = 1;
        noiseParams.spectrumShape = 'lowPassCornerFreq';
        noiseParams.cornerFrequencyCPD = 5;
    
        [highFrequencyStimulus, noiseNorm] = generateStimulusImage(stimulusWidthArcMin, stimulusPixelSizeArcMin,  gratingParams, noiseParams, noiseNorm, nTrials, figNo);
        highFreqFEM = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, highFrequencyStimulus);
        highFreqStabilized = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin*0, highFrequencyStimulus);
        save('highFrequency.mat', 'highFrequencyStimulus', 'highFreqFEM', 'highFreqStabilized', '-v7.3');
        
        figNo = figNo + 1;
        analyzeSpectra(highFreqFEM, highFreqStabilized, figNo);
    end
    
end





