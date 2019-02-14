function [timeAxis, photoCurrents, noisyPhotoCurrentsInstance, timeAxisConeExcitations, coneExcitations, noisyConeExcitationInstance, photocurrentSNR, coneExcitationSNR, legends] = ...
        computeNoisyReponseSNRs(adaptationPhotonRates, pulseWeberContrasts, ...
        pulseDurationSeconds, simulationTimeStepSeconds, eccentricity, noisyInstancesNum, timeWindowForSNRanalysis)
    
    % Compute the model responses for all conditions
    [~, ~, modelResponses, legends] = computeStepReponses(adaptationPhotonRates, ...
        pulseWeberContrasts, pulseDurationSeconds, simulationTimeStepSeconds, eccentricity, noisyInstancesNum);

    % Extract the noisy photocurrents
    nAdaptationLevels = size(modelResponses,1);
    mPulseStrengths = size(modelResponses,2);
    
    % Initialize
    photocurrentSNR = zeros(nAdaptationLevels,mPulseStrengths);
    coneExcitationSNR = zeros(nAdaptationLevels,mPulseStrengths);
    
    for adaptationIndex = 1:nAdaptationLevels 
        for pulseStrengthIndex = 1:mPulseStrengths
            modelResponse = modelResponses{adaptationIndex, pulseStrengthIndex};
            timeAxis = modelResponse.timeAxis;
            timeAxisConeExcitations = timeAxis;
            
            % Allocate memory
            if (adaptationIndex*pulseStrengthIndex == 1)
                photoCurrents = zeros(nAdaptationLevels, mPulseStrengths, numel(timeAxis));
                noisyPhotoCurrentsInstance = zeros(nAdaptationLevels, mPulseStrengths, numel(timeAxis));
            end
            
            % Extract mean photocurrent responses
            photoCurrents(adaptationIndex, pulseStrengthIndex,:) = modelResponse.membraneCurrent;
            
            % Save the first noisy instance
            noisyPhotoCurrentsInstance(adaptationIndex, pulseStrengthIndex,:) = modelResponse.noisyMembraneCurrents(1,:);
           
            % Compute the SNR based on the noise traces (signal-mean signal) and the mean response modulation
            noiseTraces = bsxfun(@minus, modelResponse.noisyMembraneCurrents, modelResponse.membraneCurrent);
            responseModulation = modelResponse.membraneCurrent - modelResponse.membraneCurrent(1);
            
            % Compute the SNR for the photocurrents
            displayVectors = ~true;
            photocurrentSNR(adaptationIndex, pulseStrengthIndex) = ...
                computeSNR(responseModulation, noiseTraces, displayVectors, timeAxis, timeWindowForSNRanalysis);
                     
            
            % Do SNR computation for the cone excitation signal
            % To compute Poisson noise we need to convert photon rate to photon count. What integation time to use?
            % Do this for pulseDurationSeconds
            photonCountingIntegrationTime = pulseDurationSeconds;
            meanConeExcitationCountSignal = round(modelResponse.pRate * photonCountingIntegrationTime);

            % Obtain noisy cone excitation response instances
            noisyConeExcitationCountSignals = coneMosaic.photonNoise(repmat(meanConeExcitationCountSignal, [noisyInstancesNum 1]));
            
            % Back to photon rates
            meanConeExcitationRateSignal = meanConeExcitationCountSignal / photonCountingIntegrationTime;
            noisyConeExcitationRateSignals = noisyConeExcitationCountSignals / photonCountingIntegrationTime;
            
            % Allocate memory
            if (adaptationIndex*pulseStrengthIndex == 1)
                coneExcitations = zeros(nAdaptationLevels, mPulseStrengths, numel(timeAxis));
                noisyConeExcitationInstance = zeros(nAdaptationLevels, mPulseStrengths, numel(timeAxis));
            end
            
            % Extract mean cone excitation responses
            coneExcitations(adaptationIndex, pulseStrengthIndex,:) = meanConeExcitationCountSignal;
            
            % Save the first noisy instance
            noisyConeExcitationInstance(adaptationIndex, pulseStrengthIndex,:) = noisyConeExcitationRateSignals(1,:);
           
            % Compute the SNR based on the noise traces (signal-mean signal) and the mean response modulation
            noiseTraces = bsxfun(@minus, noisyConeExcitationRateSignals, meanConeExcitationRateSignal);
            responseModulation = meanConeExcitationRateSignal - meanConeExcitationRateSignal(1);
            
            % Compute the SNR
            coneExcitationSNR(adaptationIndex, pulseStrengthIndex) = ...
                computeSNR(responseModulation, noiseTraces, displayVectors, timeAxis, timeWindowForSNRanalysis); 
        end
    end
end

function theSNR = computeSNR(meanSignalModulation, noiseTraces, displayVectors, timeAxis, timeWindowForSNRanalysis )

    % Concatenate signals along all their instances
    noisyInstancesNum = size(noiseTraces,1);
    timeBinsForSNR = find(timeAxis>=timeWindowForSNRanalysis(1) & timeAxis <=timeWindowForSNRanalysis(2));
    
    % Only use the signals withing the timeWindowForSNRanalysis
    meanSignalModulation = meanSignalModulation(timeBinsForSNR);
    noiseTraces = noiseTraces(:,timeBinsForSNR);
    
    % Concatenate all instances along a long vector
    signal = repmat(meanSignalModulation, [1 noisyInstancesNum]);
    noise  = reshape(noiseTraces', [1 numel(noiseTraces)]); 

    % Compute the SNR on these long vectors
    theSNR = snr(signal, noise);

    if (displayVectors)
        % Just display 3 instances of the SNR signals
        displayedTimeBins = 1:(3*numel(timeBinsForSNR));
        figure(1234); clf;
        stairs(signal(displayedTimeBins), 'r-', 'LineWidth', 2); hold on;
        stairs(noise(displayedTimeBins), 'k-');
        title(sprintf('SNR = %2.2f', theSNR));
        
        % Just display the SNR mean signal  
        figure(1235); clf;
        stairs(timeAxis(timeBinsForSNR)*1000, meanSignalModulation, 'r-', 'LineWidth', 2); hold on;
    end
    
end

