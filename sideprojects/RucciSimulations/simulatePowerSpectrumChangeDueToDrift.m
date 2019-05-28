function simulatePowerSpectrumChangeDueToDrift

    % Experiment params
    
    nTrials = 8;
    emDurationSeconds = 0.2;
    
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
    
    
    % Noise params for noise only stimulus
    gratingParams.contrast = 0.0;
    noiseParams.spectrumShape = '1overF';
    
    [oneOverFnoiseStimulus, noiseNorm] = generateStimulusImage(stimulusWidthArcMin, stimulusPixelSizeArcMin,  gratingParams, noiseParams, noiseNorm, nTrials, figNo);
    oneOverFnoiseFEM = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, oneOverFnoiseStimulus);
    oneOverFnoiseStabilized = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin*0, oneOverFnoiseStimulus);
    figNo = figNo + 1;
    analyzeSpectra(oneOverFnoiseFEM, oneOverFnoiseStabilized, figNo);
    pause;
    
    
    % Generate low frequency stimulus with fixational EM
    % Params for low frequency stimulus
    gratingParams.sfCPD = 4;
    gratingParams.contrast = 1;
    
    noiseParams.spectrumShape = 'highPassCornerFreq';
    noiseParams.cornerFrequencyCPD = 10;
    
    [lowFrequencyStimulus, noiseNorm] = generateStimulusImage(stimulusWidthArcMin, stimulusPixelSizeArcMin,  gratingParams, noiseParams, noiseNorm, nTrials, figNo);
    lowFreqFEM = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, lowFrequencyStimulus);
    lowFreqStabilized = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin*0, lowFrequencyStimulus);
    figNo = figNo + 1;
    analyzeSpectra(lowFreqFEM, lowFreqStabilized, figNo);
    
    
    
    
    % Generate high frequency stimulus without fixational EM
    
    % Noise params for highfrequency stimulus
    gratingParams.sfCPD = 11;
    gratingParams.contrast = 1;
    noiseParams.spectrumShape = 'lowPassCornerFreq';
    noiseParams.cornerFrequencyCPD = 5;
    
%     figNo = figNo + 1;
%     highFrequencyStimulus = generateStimulusImage(stimulusWidthArcMin, stimulusPixelSizeArcMin,  gratingParams, noiseParams, noiseNorm, nTrials, figNo);
%     highFreqFEM = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin*0.5, highFrequencyStimulus);
%     highFreqStabilized = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin*0, highFrequencyStimulus);
%     figNo = figNo + 1;
%     analyzeSpectra(highFreqFEM, highFreqStabilized, figNo);
    pause
    
    

    
    % Noise params for noise only stimulus
    %gratingParams.contrast = 0.0;
    %noiseParams.spectrumShape = '1overF';
    
    % Generate stimulus
    %figNo = 3;
    %noiseOnlyStimulus = generateStimulusImage(stimulusWidthArcMin, stimulusPixelSizeArcMin,  gratingParams, noiseParams, nTrials, figNo);
    
 
    % Noise params for grating only stimulus
    %gratingParams.contrast = 1.0;
    %noiseParams.spectrumShape = 'none';
    
    % Generate stimulus
    %figNo = 4;
    %gratingOnlyStimulus = generateStimulusImage(stimulusWidthArcMin, stimulusPixelSizeArcMin,  gratingParams, noiseParams, nTrials, figNo);
    
    
    
end





