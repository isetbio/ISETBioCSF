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
    gratingParams.sfCPD = [];
    
    % Noise params
    noiseParams.spectrumShape = '';
    noiseParams.cornerFrequencyCPD = [];
    noiseParams.steepness = 60;
    
    % Params for low frequency stimulus
    gratingParams.sfCPD = 4;
    noiseParams.spectrumShape = 'highPassCornerFreq';
    noiseParams.cornerFrequencyCPD = 10;
    
    % Generate stimulus
    figNo = 1;
    lowFrequencyStimulus = generateStimulusImage(stimulusWidthArcMin, stimulusPixelSizeArcMin,  gratingParams, noiseParams, nTrials, figNo);
    

    % Noise params for highfrequency stimulus
    gratingParams.sfCPD = 11;
    noiseParams.spectrumShape = 'lowPassCornerFreq';
    noiseParams.cornerFrequencyCPD = 5;
    
    % Generate stimulus
    figNo = 2;
    highFrequencyStimulus = generateStimulusImage(stimulusWidthArcMin, stimulusPixelSizeArcMin,  gratingParams, noiseParams, nTrials, figNo);
    
    
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
    
    
    generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, highFrequencyStimulus, 100);
    generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin*0, highFrequencyStimulus, 200);
    
    generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, lowFrequencyStimulus, 300);
    generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin*0, lowFrequencyStimulus, 400);
    
    
    
end





