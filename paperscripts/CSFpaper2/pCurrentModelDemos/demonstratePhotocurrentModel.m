function demonstratePhotocurrentModel
    
    % eccentricity for biophysical model
    eccentricity = 'foveal';
    
    % time step for biophysical model
    simulationTimeStepSeconds = 0.1/1000;
   
    
    doImpulseResponseAnalysis = ~true;
    doOnOffAsymmetryAnalysis = ~true;
    doSignalToNoiseAnalysis = true;
    doStepResponseAnalysis = ~true;                     % used to illustrate the different model components
    doBipolarStepResponseAnalysis = ~true;
    
    figNo = 0;
    
    % Compute and visualize the pCurrent impulse responses at different adaptation levels
    if (doImpulseResponseAnalysis)
        figNo = demoImpulseResponses(eccentricity, simulationTimeStepSeconds, figNo);
    end
    
    % Compute and visualze the pCurrent step responses at differentadaptation levels
    if (doStepResponseAnalysis)
        adaptationPhotonRate = 30 / 5 *  1000; % corresponding to typical 5 ms count of 30 cone excitations
        pulseWeberContrast = 0.90;
        figNo = demoStepResponses(eccentricity, simulationTimeStepSeconds, adaptationPhotonRate, pulseWeberContrast, figNo);
    end
    % Compute and visualze the pCurrent step responses at differentadaptation levels
    if (doBipolarStepResponseAnalysis)
        figNo = demoBipolarStepResponses(eccentricity, simulationTimeStepSeconds, figNo);
    end
    
    adaptationPhotonRates = [50 100 200 400 600 1000 1800 3000 5000 10000 18000];
    spontaneousIsomerizationRate = 50;
    
    noisyInstancesNum = 10;
        
    if (doOnOffAsymmetryAnalysis)
        figNo = demoOnOFFResponses(eccentricity, simulationTimeStepSeconds, spontaneousIsomerizationRate, ...
            adaptationPhotonRates, pulseDurationSeconds, figNo);
    end

    if (doSignalToNoiseAnalysis)
        
        pulseDurationSeconds = 50/1000; integrationTimeSeconds = 300/1000;
        figNo = demoSignalToNoiseRatio(eccentricity, simulationTimeStepSeconds, spontaneousIsomerizationRate, ...
            adaptationPhotonRates, pulseDurationSeconds, integrationTimeSeconds,  noisyInstancesNum, figNo);
        
        pulseDurationSeconds = 50/1000; integrationTimeSeconds = 200/1000;
        figNo = demoSignalToNoiseRatio(eccentricity, simulationTimeStepSeconds, spontaneousIsomerizationRate, ...
            adaptationPhotonRates, pulseDurationSeconds, integrationTimeSeconds,  noisyInstancesNum, figNo);
        
        pulseDurationSeconds = 50/1000; integrationTimeSeconds = 100/1000;
        figNo = demoSignalToNoiseRatio(eccentricity, simulationTimeStepSeconds, spontaneousIsomerizationRate, ...
            adaptationPhotonRates, pulseDurationSeconds, integrationTimeSeconds,  noisyInstancesNum, figNo);
        
        pulseDurationSeconds = 50/1000; integrationTimeSeconds = 50/1000;
        figNo = demoSignalToNoiseRatio(eccentricity, simulationTimeStepSeconds, spontaneousIsomerizationRate, ...
            adaptationPhotonRates, pulseDurationSeconds, integrationTimeSeconds,  noisyInstancesNum, figNo);
        
        pulseDurationSeconds = 50/1000; integrationTimeSeconds = 20/1000;
        figNo = demoSignalToNoiseRatio(eccentricity, simulationTimeStepSeconds, spontaneousIsomerizationRate, ...
            adaptationPhotonRates, pulseDurationSeconds, integrationTimeSeconds,  noisyInstancesNum, figNo);
        
    end
    
end

function figNo = demoOnOFFResponses(eccentricity, simulationTimeStepSeconds, spontaneousIsomerizationRate, adaptationPhotonRates, pulseDurationSeconds, figNo)
    % compute the responses the for different step Weber contrasts
    pulseWeberContrasts = logspace(log10(0.05), log10(0.95), 10);  % [0.15 0.25 0.5 0.75 1.0];
    interpulseIntervalSeconds = 600/1000;
    
    % Define the stim params struct
    constantStimParams = struct(...
        'type', 'on_off_pulses', ...                            % type of stimulus
        'pulseDurationSeconds', pulseDurationSeconds, ...       % pulse duration in seconds
        'interpulseIntervalSeconds', interpulseIntervalSeconds, ...                     % interpular interval in seconds
        'totalDurationSeconds', 1.3, ...                          % total duration of the stimulus
        'timeSampleSeconds', simulationTimeStepSeconds ...
    );

    % no noisy  instances
    instancesNum = 0;
    %false to visualize the internal model components
    useDefaultPhotocurrentImplementation = true;
    
    [timeAxis, onOffResponses, modelResponses, legends] = ...
        computeStepReponses(constantStimParams, adaptationPhotonRates, pulseWeberContrasts, ...
        pulseDurationSeconds, spontaneousIsomerizationRate, eccentricity, instancesNum, useDefaultPhotocurrentImplementation);
    
    % Plot the different On/OFF responses
    figNo = figNo + 1;
    plotOnOffResponses(timeAxis, onOffResponses, modelResponses, adaptationPhotonRates, pulseWeberContrasts, interpulseIntervalSeconds+pulseDurationSeconds, figNo);
end

function figNo = demoSignalToNoiseRatio(eccentricity, simulationTimeStepSeconds, spontaneousIsomerizationRate, adaptationPhotonRates, pulseDurationSeconds, integrationTimeSeconds, noisyInstancesNum, figNo)

    adaptationPhotonRates = round(adaptationPhotonRates/10)*10;
    idx = find(adaptationPhotonRates>100);
    adaptationPhotonRates(idx) = round(adaptationPhotonRates(idx)/100)*100;
    pulseWeberContrasts   = [0.025 0.05 0.1 0.2 0.4 0.8]; % logspace(log10(0.03), log10(0.30), nContrastLevels);
    
    
    % analyze SNR using the [0 200] time period
    timeWindowForSNRanalysis =  [0 200]/1000;
    
    % analyze SNR using the peak response
    %timeWindowForSNRanalysis = [];
    
    % analyze SNR using a +/- 0.5*pulseDurationSeconds around the peak time
    timeWindowForSNRanalysis = integrationTimeSeconds; %pulseDurationSeconds;
    
    displaySNRVectors = true;
    transformDecibelsToRatios = true;
    
    [timeAxis, photoCurrents, noisyPhotoCurrentsInstance, ...
     timeAxisConeExcitations, coneExcitationRates, noisyConeExcitationRateInstance, ...
     photocurrentSNR, coneExcitationSNR, legends] = computeNoisyReponseSNRs(...
        adaptationPhotonRates, pulseWeberContrasts, ...
        pulseDurationSeconds, simulationTimeStepSeconds, spontaneousIsomerizationRate, eccentricity, noisyInstancesNum, timeWindowForSNRanalysis, displaySNRVectors);
    
    
    % Plot the different step responses
    figNo = figNo + 1;
    SNRLims = [-35 30]; % [-50 35];
    SNRTicks = -30:5:50;
    
    plotSignalToNoiseResults(timeAxis, photoCurrents, noisyPhotoCurrentsInstance, ...
        timeAxisConeExcitations, coneExcitationRates, noisyConeExcitationRateInstance, ...
        photocurrentSNR, coneExcitationSNR, transformDecibelsToRatios, adaptationPhotonRates, ...
        pulseWeberContrasts, SNRLims, SNRTicks, spontaneousIsomerizationRate, pulseDurationSeconds, integrationTimeSeconds, legends, figNo); 
end


function figNo = demoImpulseResponses(eccentricity, simulationTimeStepSeconds, figNo)
    % Examined adaptation levels (photons/cone/sec)
    adaptationPhotonRates = [60 600 6000]; % [0 300 1000 3000 10000];
    impulseDurationSeconds = simulationTimeStepSeconds;
    photonCountDuringImpulse = 1;
    %false to visualize the internal model components
    useDefaultPhotocurrentImplementation = true;
    
    % Compute the impulse response at different adaptation levels
    [timeAxis, impulseResponses, temporalFrequencyAxis, impulseResponseSpectra, modelResponses, legends] = ...
        computeImpulseReponses(impulseDurationSeconds, photonCountDuringImpulse, adaptationPhotonRates, ...
        simulationTimeStepSeconds, eccentricity, useDefaultPhotocurrentImplementation);
    
    % Plot the impulse response at different adaptation levels
    figNo = figNo + 1;
    plotFrequencySpectra = ~true;
    plotImpulseResponses(timeAxis, impulseResponses, temporalFrequencyAxis, impulseResponseSpectra, adaptationPhotonRates, plotFrequencySpectra, legends, figNo);
    
    % Plot the model response to the impulse stimuli
    if (~useDefaultPhotocurrentImplementation)
        figNo = figNo + 1;
        plotModelResponses(modelResponses, legends,figNo);
    end
end


function figNo = demoStepResponses(eccentricity, simulationTimeStepSeconds, adaptationPhotonRate, pulseWeberContrast,figNo)
    % Examined adaptation levels (photons/cone/sec)
    adaptationPhotonRates = adaptationPhotonRate;
    stepDurationSeconds = 500/1000;
    photonCountDuringImpulse = adaptationPhotonRates*stepDurationSeconds*pulseWeberContrast;
    
    %false to visualize the internal model components
    useDefaultPhotocurrentImplementation = true;
    
    % Compute the impulse response at different adaptation levels
    [timeAxis, impulseResponses, temporalFrequencyAxis, impulseResponseSpectra, modelResponses, legends] = ...
        computeImpulseReponses(stepDurationSeconds, photonCountDuringImpulse, adaptationPhotonRates, simulationTimeStepSeconds, eccentricity, useDefaultPhotocurrentImplementation);
    
    % Plot the model response to the impulse stimuli
    figNo = figNo + 1;
    plotModelResponses(modelResponses, legends, []);
end


function figNo = demoBipolarStepResponses(eccentricity, simulationTimeStepSeconds, figNo)

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
    %false to visualize the internal model components
    useDefaultPhotocurrentImplementation = true;
    
    [timeAxis, stepResponses, modelResponses, legends] = ...
        computeStepReponses(constantStimParams, adaptationPhotonRates, pulseWeberContrasts, ...
        pulseDurationSeconds, spontaneousIsomerizationRate, eccentricity, instancesNum, useDefaultPhotocurrentImplementation);
    
    % Plot the different step responses
    figNo = figNo + 1;
    plotStepResponses(timeAxis, stepResponses, adaptationPhotonRates, pulseWeberContrasts, legends, figNo);
end