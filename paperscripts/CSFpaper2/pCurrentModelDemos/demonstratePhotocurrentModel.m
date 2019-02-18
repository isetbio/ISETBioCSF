function demonstratePhotocurrentModel
    figNo = 0;
    eccentricity = 'foveal';
    simulationTimeStepSeconds = 0.05/1000;
    
    doImpulseResponseAnalysis = ~true;
    doStepResponseAnalysis = ~true;
    doSignalToNoiseAnalysis = true;
    
    % Compute and visualize the pCurrent impulse responses at different adaptation levels
    if (doImpulseResponseAnalysis)
        figNo = demoImpulseResponses(eccentricity, simulationTimeStepSeconds, figNo);
    end
    
    % Compute and visualze the pCurrent step responses at differentadaptation levels
    if (doStepResponseAnalysis)
        figNo = demoStepResponses(eccentricity, simulationTimeStepSeconds, figNo);
    end
    
    if (doSignalToNoiseAnalysis)
        figNo = demoSignalToNoiseRatio(eccentricity, simulationTimeStepSeconds, figNo);
    end
    
end

function figNo = demoSignalToNoiseRatio(eccentricity, simulationTimeStepSeconds, figNo)

    % Examined adaptation levels (photons/cone/sec)
    nAdaptationLevels = 5;
    nContrastLevels = 7;
    adaptationPhotonRates = logspace(log10(0.1), log10(10), nAdaptationLevels) * 1000;
    pulseWeberContrasts   = logspace(log10(0.02),log10(0.90),nContrastLevels);
    pulseDurationSeconds = 100/1000;
    
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
        pulseDurationSeconds, simulationTimeStepSeconds, eccentricity, noisyInstancesNum, timeWindowForSNRanalysis, displaySNRVectors);
    
    
    % Plot the different step responses
    figNo = figNo + 1;
    plotSignalToNoiseResults(timeAxis, photoCurrents, noisyPhotoCurrentsInstance, ...
        timeAxisConeExcitations, coneExcitationRates, noisyConeExcitationRateInstance, ...
        photocurrentSNR, coneExcitationSNR, transformDecibelsToRatios, adaptationPhotonRates, pulseWeberContrasts, legends, figNo); 
end


function figNo = demoImpulseResponses(eccentricity, simulationTimeStepSeconds, figNo)
    % Examined adaptation levels (photons/cone/sec)
    adaptationPhotonRates = [0 300 1000 3000 10000];
    
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
    
    [timeAxis, stepResponses, modelResponses, legends] = ...
        computeStepReponses(adaptationPhotonRates, pulseWeberContrasts, ...
        pulseDurationSeconds, simulationTimeStepSeconds, eccentricity);
    
    % Plot the different step responses
    figNo = figNo + 1;
    plotStepResponses(timeAxis, stepResponses, adaptationPhotonRates, pulseWeberContrasts, legends, figNo);
end