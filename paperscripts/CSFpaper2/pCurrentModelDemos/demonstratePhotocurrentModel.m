function demonstratePhotocurrentModel
    
    % eccentricity for biophysical model
    eccentricity = 'foveal';
    
    % time step for biophysical model
    simulationTimeStepSeconds = 0.1/1000;
    
    % Spontaneous photoisomerization rate (100 R*/c/sec)
    spontaneousIsomerizationRate = 100;
    
    doImpulseResponseAnalysis = true;
    doOnOffAsymmetryAnalysis = ~true;
    doStepResponseAnalysis = ~true;
    doSignalToNoiseAnalysis = true;
    
    
    figNo = 0;
    
    % Compute and visualize the pCurrent impulse responses at different adaptation levels
    if (doImpulseResponseAnalysis)
        figNo = demoImpulseResponses(eccentricity, simulationTimeStepSeconds, figNo);
    end
    
    if (doOnOffAsymmetryAnalysis)
        figNo = demoOnOFFResponses(eccentricity, simulationTimeStepSeconds, figNo);
    end
    
    % Compute and visualze the pCurrent step responses at differentadaptation levels
    if (doStepResponseAnalysis)
        figNo = demoStepResponses(eccentricity, simulationTimeStepSeconds, figNo);
    end
    
    if (doSignalToNoiseAnalysis)
        figNo = demoSignalToNoiseRatio(eccentricity, simulationTimeStepSeconds, spontaneousIsomerizationRate, figNo);
    end
    
end

function figNo = demoOnOFFResponses(eccentricity, simulationTimeStepSeconds, figNo)
    % Compute the step responses at different adaptation levels
    adaptationPhotonRates = [1000 10000 13000];
    pulseDurationSeconds = 1000/1000;
    pulseWeberContrasts = [0.12 0.24 0.48 0.96];
   
    
    % Define the stim params struct
    constantStimParams = struct(...
        'type', 'on_off_pulses', ...                            % type of stimulus
        'pulseDurationSeconds', pulseDurationSeconds, ...       % pulse duration in seconds
        'interpulseIntervalSeconds', 1, ...                     % interpular interval in seconds
        'totalDurationSeconds', 4, ...                          % total duration of the stimulus
        'timeSampleSeconds', simulationTimeStepSeconds ...
    );

    spontaneousIsomerizationRate = 0;
    instancesNum = 1;
    
    [timeAxis, onOffResponses, modelResponses, legends] = ...
        computeStepReponses(constantStimParams, adaptationPhotonRates, pulseWeberContrasts, ...
        pulseDurationSeconds, spontaneousIsomerizationRate, eccentricity, instancesNum);
    
    % Plot the different On/OFF responses
    figNo = figNo + 1;
    plotOnOffResponses(timeAxis, onOffResponses, modelResponses, adaptationPhotonRates, pulseWeberContrasts, legends, figNo);
end

function figNo = demoSignalToNoiseRatio(eccentricity, simulationTimeStepSeconds, spontaneousIsomerizationRate, figNo)

    % Examined adaptation levels (photons/cone/sec)
    nAdaptationLevels = 10;
    nContrastLevels = 3;
    adaptationPhotonRates = logspace(log10(0.03), log10(10),   nAdaptationLevels) * 1000;
    adaptationPhotonRates = round(adaptationPhotonRates/10)*10;
    idx = find(adaptationPhotonRates>100);
    adaptationPhotonRates(idx) = round(adaptationPhotonRates(idx)/100)*100;
    pulseWeberContrasts   = [0.016 0.1 0.6]; % logspace(log10(0.03), log10(0.30), nContrastLevels);
    pulseDurationSeconds  = 100/1000;
    
    noisyInstancesNum = 100;
    
    % analyze SNR using the [0 200] time period
    timeWindowForSNRanalysis =  [0 200]/1000;
    
    % analyze SNR using the peak response
    timeWindowForSNRanalysis = [];
    
    % analyze SNR using a +/- 50 msec window around the peak time
    timeWindowForSNRanalysis = 100/1000;
    
    displaySNRVectors = true;
    transformDecibelsToRatios = true;
    
    [timeAxis, photoCurrents, noisyPhotoCurrentsInstance, ...
     timeAxisConeExcitations, coneExcitationRates, noisyConeExcitationRateInstance, ...
     photocurrentSNR, coneExcitationSNR, legends] = computeNoisyReponseSNRs(...
        adaptationPhotonRates, pulseWeberContrasts, ...
        pulseDurationSeconds, simulationTimeStepSeconds, spontaneousIsomerizationRate, eccentricity, noisyInstancesNum, timeWindowForSNRanalysis, displaySNRVectors);
    
    
    % Plot the different step responses
    figNo = figNo + 1;
    plotSignalToNoiseResults(timeAxis, photoCurrents, noisyPhotoCurrentsInstance, ...
        timeAxisConeExcitations, coneExcitationRates, noisyConeExcitationRateInstance, ...
        photocurrentSNR, coneExcitationSNR, transformDecibelsToRatios, adaptationPhotonRates, pulseWeberContrasts, legends, figNo); 
end


function figNo = demoImpulseResponses(eccentricity, simulationTimeStepSeconds, figNo)
    % Examined adaptation levels (photons/cone/sec)
    adaptationPhotonRates = [0]; % [0 300 1000 3000 10000];
    
    % Compute the impulse response at different adaptation levels
    [timeAxis, impulseResponses, temporalFrequencyAxis, impulseResponseSpectra, modelResponses, legends] = ...
        computeImpulseReponses(adaptationPhotonRates, simulationTimeStepSeconds, eccentricity);
    
    % Plot the impulse response at different adaptation levels
    figNo = figNo + 1;
    plotImpulseResponses(timeAxis, impulseResponses, temporalFrequencyAxis, impulseResponseSpectra, adaptationPhotonRates, legends, figNo);
    
    % Plot the model response to the impulse stimuli
    figNo = figNo + 1;
    plotModelResponses(modelResponses, legends,figNo);
end


function figNo = demoStepResponses(eccentricity, simulationTimeStepSeconds, figNo)

    % Compute the step responses at different adaptation levels
    adaptationPhotonRates = [100 300 1000 3000 10000];
    pulseDurationSeconds = 100/1000;
    pulseWeberContrasts = [0.125 0.25 0.5 0.75 1.0];
    
    % Define the stim params struct
    constantStimParams = struct(...
        'type', 'pulse', ...                                    % type of stimulus
        'pulseDurationSeconds', pulseDurationSeconds, ...       % pulse duration in seconds
        'totalDurationSeconds', 0.6, ...                        % total duration of the stimulus
        'timeSampleSeconds', simulationTimeStepSeconds ...
    );

    spontaneousIsomerizationRate = 0;
    instancesNum = 1;
    
    [timeAxis, stepResponses, modelResponses, legends] = ...
        computeStepReponses(constantStimParams, adaptationPhotonRates, pulseWeberContrasts, ...
        pulseDurationSeconds, spontaneousIsomerizationRate, eccentricity, instancesNum);
    
    % Plot the different step responses
    figNo = figNo + 1;
    plotStepResponses(timeAxis, stepResponses, adaptationPhotonRates, pulseWeberContrasts, legends, figNo);
end