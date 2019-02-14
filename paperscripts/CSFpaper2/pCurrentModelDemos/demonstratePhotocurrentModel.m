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
    adaptationPhotonRates = [100 300 1000 3000 10000];
    pulseDurationSeconds = 100/1000;
    pulseWeberContrasts = logspace(log10(0.01),log10(1.0),5);
    noisyInstancesNum = 1000;
    
    [timeAxis, photoCurrents, noisyPhotoCurrentsInstance, ...
     timeAxisConeExcitations, coneExcitations, noisyConeExcitationInstance, photocurrentSNR, coneExcitationSNR, legends] = ...
        computeNoisyPhotocurrentReponses(adaptationPhotonRates, pulseWeberContrasts, ...
        pulseDurationSeconds, simulationTimeStepSeconds, eccentricity, noisyInstancesNum);
    
    
    % Plot the different step responses
    figNo = figNo + 1;
    plotSignalToNoiseResults(timeAxis, photoCurrents, noisyPhotoCurrentsInstance, ...
        timeAxisConeExcitations, coneExcitations, noisyConeExcitationInstance, ...
        photocurrentSNR, coneExcitationSNR, adaptationPhotonRates, pulseWeberContrasts, legends, figNo); 
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