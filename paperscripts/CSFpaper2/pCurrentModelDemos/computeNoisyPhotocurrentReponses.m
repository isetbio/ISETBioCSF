function [timeAxis, photoCurrents, noisyPhotoCurrentsInstance, timeAxisConeExcitations, coneExcitations, noisyConeExcitationInstance, photocurrentSNR, coneExcitationSNR, legends] = ...
        computeNoisyPhotocurrentReponses(adaptationPhotonRates, pulseWeberContrasts, ...
        pulseDurationSeconds, simulationTimeStepSeconds, eccentricity, noisyInstancesNum)
    
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
            
            % Allocate memory
            if (adaptationIndex*pulseStrengthIndex == 1)
                photoCurrents = zeros(nAdaptationLevels, mPulseStrengths, numel(timeAxis));
                noisyPhotoCurrentsInstance = zeros(nAdaptationLevels, mPulseStrengths, numel(timeAxis));
            end
            
            % Extract mean photocurrent responses
            photoCurrents(adaptationIndex, pulseStrengthIndex,:) = modelResponse.membraneCurrent;
            
            % Save the first noisy instance
            noisyPhotoCurrentsInstance(adaptationIndex, pulseStrengthIndex,:) = modelResponse.noisyMembraneCurrents(1,:);
           
            % Compute the SNR for each noisy instance
            noisyTraces = bsxfun(@minus, modelResponse.noisyMembraneCurrents, modelResponse.membraneCurrent);
            theSNRs = zeros(1, noisyInstancesNum);
            for instance = 1:noisyInstancesNum
                theSNRs(instance) = snr(modelResponse.membraneCurrent, noisyTraces(instance,:));
            end
            
            % Compute mean SNR across all instances
            photocurrentSNR(adaptationIndex, pulseStrengthIndex) = mean(theSNRs);
            
            % Do similar SNR computation for the cone excitation signal
            
            % Rates to counts. What integation time to use?
            meanConeExcitationCountSignal = modelResponse.pRate * pulseDurationSeconds;
            
            % subsample?
            dt = 5/1000; % assume 5 msec binwidth
            stepInc = round(dt/(timeAxis(2)-timeAxis(1)));
            idx = 1:stepInc:numel(timeAxis);
            timeAxisConeExcitations = timeAxis(idx);
            meanConeExcitationCountSignal = meanConeExcitationCountSignal(idx);
            noisyConeExcitationCountSignals = coneMosaic.photonNoise(repmat(meanConeExcitationCountSignal, [noisyInstancesNum 1]));
            
            %back to photon rates
            meanConeExcitationCountSignal = meanConeExcitationCountSignal / pulseDurationSeconds;
            noisyConeExcitationCountSignals = noisyConeExcitationCountSignals / pulseDurationSeconds;
            
            % Allocate memory
            if (adaptationIndex*pulseStrengthIndex == 1)
                coneExcitations = zeros(nAdaptationLevels, mPulseStrengths, numel(timeAxisConeExcitations));
                noisyConeExcitationInstance = zeros(nAdaptationLevels, mPulseStrengths, numel(timeAxisConeExcitations));
            end
            
            % Extract mean cone excitation responses
            coneExcitations(adaptationIndex, pulseStrengthIndex,:) = meanConeExcitationCountSignal;
            
            % Save the first noisy instance
            noisyConeExcitationInstance(adaptationIndex, pulseStrengthIndex,:) = noisyConeExcitationCountSignals(1,:);
           
            % Compute the SNR for each noisy cone excitation instance
            noisyTraces = bsxfun(@minus, noisyConeExcitationCountSignals, meanConeExcitationCountSignal);
            for instance = 1:noisyInstancesNum
                theSNRs(instance) = snr(meanConeExcitationCountSignal, noisyTraces(instance,:));
            end
            
            % Compute mean SNR across all instances
            coneExcitationSNR(adaptationIndex, pulseStrengthIndex) = mean(theSNRs);
        end
    end
end

