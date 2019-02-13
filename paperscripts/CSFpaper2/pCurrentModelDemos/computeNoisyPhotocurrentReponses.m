function [timeAxis, photoCurrents, noisyPhotoCurrentsInstance, photocurrentSNR, coneExcitationSNR, legends] = ...
        computeNoisyPhotocurrentReponses(adaptationPhotonRates, pulseWeberContrasts, ...
        pulseDurationSeconds, simulationTimeStepSeconds, eccentricity, noisyInstancesNum)
    
    % Compute the model responses for all conditions
    [~, ~, modelResponses, legends] = computeStepReponses(adaptationPhotonRates, ...
        pulseWeberContrasts, pulseDurationSeconds, simulationTimeStepSeconds, eccentricity, noisyInstancesNum);

    % Extract the noisy photocurrents
    nAdaptationLevels = size(modelResponses,1);
    mPulseStrengths = size(modelResponses,2);
    
    for adaptationIndex = 1:nAdaptationLevels 
        for pulseStrengthIndex = 1:mPulseStrengths
            modelResponse = modelResponses{adaptationIndex, pulseStrengthIndex};
            timeAxis = modelResponse.timeAxis;
            if (adaptationIndex*pulseStrengthIndex == 1)
                photoCurrents = zeros(nAdaptationLevels, mPulseStrengths, numel(modelResponse.membraneCurrent));
                noisyPhotoCurrentsInstance = zeros(nAdaptationLevels, mPulseStrengths, numel(modelResponse.membraneCurrent));
            end
            photoCurrents(adaptationIndex, pulseStrengthIndex,:) = modelResponse.membraneCurrent;
            noisyPhotoCurrentsInstance(adaptationIndex, pulseStrengthIndex,:) = modelResponse.noisyMembraneCurrents(1,:);     % first instance
            
            % Standard deviation of pCurrent across all instances, all time samples
            sigmaPhotocurrent = std(modelResponse.noisyMembraneCurrents(:)); 
            
            % Peak modulation of pCurrent
            meanPhotocurrentTrace = photoCurrents(adaptationIndex, pulseStrengthIndex,:);
            meanPhotocurrentPeak = max(abs(meanPhotocurrentTrace-meanPhotocurrentTrace(1)));
            
            % SNR of pCurrent = (peak modulation / std)^2
            photocurrentSNR(adaptationIndex, pulseStrengthIndex) = (meanPhotocurrentPeak/sigmaPhotocurrent)^2;
            
            % Do similar SNR computation for the cone excitation signal
            % Find time bins corresponding to the pulse duration
            idx = find(modelResponse.pRate ~= adaptationPhotonRates(adaptationIndex));
            
            % Extract photon rate during the pulse
            pulsePhotonRate = mean(modelResponse.pRate(idx));
            
            % And corresponding photon count
            meanPhotonCountDuringPulse = pulsePhotonRate * pulseDurationSeconds;
            
            % Compute 10000 instances of Poisson noise with mean = meanPhotonCountDuringPulse
            photonCountDistributionDuringPulse = coneMosaic.photonNoise(meanPhotonCountDuringPulse*ones(1,10000));
            
            % Compute standard deviation of these samples
            sigmaPhotons = std(photonCountDistributionDuringPulse);
            
            % SNR of cone excitation = (meanPhotonCountDuringPulse / std of photons)^2
            coneExcitationSNR(adaptationIndex, pulseStrengthIndex) = (meanPhotonCountDuringPulse/sigmaPhotons)^2;
        end
    end
end

