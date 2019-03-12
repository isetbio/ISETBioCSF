function [timeAxis, photoCurrents, noisyPhotoCurrentsInstance, ...
    timeAxisConeExcitations, coneExcitationRates, noisyConeExcitationRateInstance, ...
    photocurrentSNR, coneExcitationSNR, legends] = ...
        computeNoisyReponseSNRs(adaptationPhotonRates, pulseWeberContrasts, ...
        pulseDurationSeconds, simulationTimeStepSeconds, spontaneousIsomerizationRate, eccentricity, noisyInstancesNum, ...
        timeWindowForSNRanalysis, displaySNRVectors)
    
    % Define the stim params struct
    constantStimParams = struct(...
        'type', 'pulse', ...                                    % type of stimulus
        'pulseDurationSeconds', pulseDurationSeconds, ...       % pulse duration in seconds
        'totalDurationSeconds', 0.4, ...                        % total duration of the stimulus
        'timeSampleSeconds', simulationTimeStepSeconds ...
    );

    %false to visualize the internal model components
    useDefaultPhotocurrentImplementation = true;
    
    % Compute the model responses for all conditions
    [~, ~, modelResponses, legends] = computeStepReponses(constantStimParams, adaptationPhotonRates, ...
        pulseWeberContrasts, pulseDurationSeconds, spontaneousIsomerizationRate, eccentricity, noisyInstancesNum, useDefaultPhotocurrentImplementation);

    % Extract the noisy photocurrents
    nAdaptationLevels = size(modelResponses,1);
    mPulseStrengths = size(modelResponses,2);
    
    % Initialize
    photocurrentSNR = zeros(nAdaptationLevels,mPulseStrengths);
    coneExcitationSNR = zeros(nAdaptationLevels,mPulseStrengths);
    
    for adaptationIndex = 1:nAdaptationLevels 
        for pulseStrengthIndex = 1:mPulseStrengths
            modelResponse = modelResponses{adaptationIndex, pulseStrengthIndex};
            
            % Compute cone excitation signals fro stimulus photon rate
            noisyConeExcitationRateSignals = ...
                computeConeExcitationSignalsFromStimulusPhotonRate(modelResponse.pRate, spontaneousIsomerizationRate, pulseDurationSeconds, noisyInstancesNum);
            meanConeExcitationRateSignal = modelResponse.pRate;
             
            timeAxis = modelResponse.timeAxis;
            timeAxisConeExcitations = timeAxis;
            
            % Allocate memory
            if (adaptationIndex*pulseStrengthIndex == 1)
                photoCurrents   = zeros(nAdaptationLevels, mPulseStrengths, numel(timeAxis));
                coneExcitationRates = zeros(nAdaptationLevels, mPulseStrengths, numel(timeAxis));
                noisyPhotoCurrentsInstance  = zeros(nAdaptationLevels, mPulseStrengths, numel(timeAxis));
                noisyConeExcitationRateInstance = zeros(nAdaptationLevels, mPulseStrengths, numel(timeAxis));
            end
            
            
            % Do SNR computation for the cone excitation signal
            % Extract mean cone excitation responses
            coneExcitationRates(adaptationIndex, pulseStrengthIndex,:) = meanConeExcitationRateSignal;
            
            % Save the first noisy instance
            noisyConeExcitationRateInstance(adaptationIndex, pulseStrengthIndex,:) = noisyConeExcitationRateSignals(1,:);
           
            % Compute the SNR based on the noise traces (signal-mean signal) and the mean response modulation
            if (any(noisyConeExcitationRateSignals(:)<0))
                fprintf('What the heck? ''noisyConeExcitationRateSignals'' has negative elements.\n');
                error('this shouldnt happen');
            end
            noiseTraces = bsxfun(@minus, noisyConeExcitationRateSignals, meanConeExcitationRateSignal);
            responseModulation = meanConeExcitationRateSignal - meanConeExcitationRateSignal(1);
            
            % Compute the SNR
            coneExcitationSNR(adaptationIndex, pulseStrengthIndex) = ...
                computeSNR(responseModulation, noiseTraces, displaySNRVectors, timeAxis, ...
                timeWindowForSNRanalysis, 'cone excitations'); 
            
            
            % Do SNR computation for the photocurrent signal
            % Extract mean photocurrent responses
            photoCurrents(adaptationIndex, pulseStrengthIndex,:) = modelResponse.membraneCurrent;
            
            % Save the first noisy instance
            noisyPhotoCurrentsInstance(adaptationIndex, pulseStrengthIndex,:) = modelResponse.noisyMembraneCurrents(1,:);
           
            % Compute the SNR based on the noise traces (signal-mean signal) and the mean response modulation
            noiseTraces = bsxfun(@minus, modelResponse.noisyMembraneCurrents, modelResponse.membraneCurrent);
            responseModulation = modelResponse.membraneCurrent - modelResponse.membraneCurrent(1);
            
            % Compute the SNR for the photocurrents
            photocurrentSNR(adaptationIndex, pulseStrengthIndex) = ...
                computeSNR(responseModulation, noiseTraces, displaySNRVectors, timeAxis, ...
                timeWindowForSNRanalysis, 'photocurrent');
                     
        end
    end
end


function noisyConeExcitationRateSignals = computeConeExcitationSignalsFromStimulusPhotonRate(pRate, spontaneousIsomerizationRate, pulseDurationSeconds, noisyInstancesNum)
            
    % To compute Poisson noise we need to convert photon rate to photon count. What integation time to use?
    % Do this for pulseDurationSeconds
    photonCountingIntegrationTime = pulseDurationSeconds;
    meanConeExcitationCountSignal = round((pRate+spontaneousIsomerizationRate) * photonCountingIntegrationTime);

    % Obtain noisy cone excitation response instances
    noisyConeExcitationCountSignals = coneMosaic.photonNoise(repmat(meanConeExcitationCountSignal, [noisyInstancesNum 1]));
    
    % Noisy photon rates
    noisyConeExcitationRateSignals = noisyConeExcitationCountSignals / photonCountingIntegrationTime;          
end


function theSNR = computeSNR(meanSignalModulation, noiseTraces, displaySNRVectors, timeAxis, timeWindowForSNRanalysis, signalName)

    % Concatenate signals along all their instances
    noisyInstancesNum = size(noiseTraces,1);
    if (numel(timeWindowForSNRanalysis) == 2)
        timeBinsForSNR = find(timeAxis>=timeWindowForSNRanalysis(1) & timeAxis <=timeWindowForSNRanalysis(2));
        duration = timeAxis(timeBinsForSNR(end))-timeAxis(timeBinsForSNR(1));
        fprintf('SNR for %s will be performed over t = [%2.2f - %2.2f] (duration:%2.2f) msec\n', signalName, timeAxis(timeBinsForSNR(1))*1000, timeAxis(timeBinsForSNR(end))*1000, duration*1000);
    else
        [~,timeBinOfPeakSignalModulation] = max(abs(meanSignalModulation));
        if (isempty(timeWindowForSNRanalysis))
            timeBinsForSNR = timeBinOfPeakSignalModulation;
            fprintf('SNR for %s will be performed at t = %2.2f  msec (peak)\n', signalName, timeAxis(timeBinsForSNR)*1000);
        else
            if (strcmp(signalName, 'cone excitations'))
                % Cone excitation windowing is simple
                t1 = timeAxis(timeBinOfPeakSignalModulation);
                t2 = t1 + timeWindowForSNRanalysis;
            else
                % Photocurrent windowing
                t1 = timeAxis(timeBinOfPeakSignalModulation) - timeWindowForSNRanalysis/2;
                t2 = timeAxis(timeBinOfPeakSignalModulation) + timeWindowForSNRanalysis/2;
                % Search around the peak to find a time point over which
                % the signal has dropped the same amount from its peak
                timeBins = find(timeAxis>=t1 & timeAxis <=t2);
                dBins = round(20/(1000*(timeAxis(2)-timeAxis(1))));
                deltaR = zeros(1,2*dBins+1);
                deltaBins = -dBins:dBins;
                
                for k = 1:numel(deltaBins)
                    if ((timeBins+deltaBins(k)>=1) && (timeBins+deltaBins(k) < numel(meanSignalModulation)))
                        r = abs(meanSignalModulation(timeBins+deltaBins(k)));
                        deltaR(k) = abs(r(1)-r(end));
                    else
                        deltaR(k) = inf;
                    end
                end
                
                [~,idx] = min(deltaR);
                % adjusted timeBinOfPeakSignalModulation
                timeBinOfPeakSignalModulation = timeBinOfPeakSignalModulation+deltaBins(idx);
                t1 = timeAxis(timeBinOfPeakSignalModulation) - timeWindowForSNRanalysis/2;
                t2 = timeAxis(timeBinOfPeakSignalModulation) + timeWindowForSNRanalysis/2;
            end
            timeBinsForSNR = find(timeAxis>=t1 & timeAxis <=t2);
            duration = timeAxis(timeBinsForSNR(end))-timeAxis(timeBinsForSNR(1));
            fprintf('SNR for %s will be performed over t = [%2.2f - %2.2f] (duration:%2.2f) msec\n', signalName, timeAxis(timeBinsForSNR(1))*1000, timeAxis(timeBinsForSNR(end))*1000, duration*1000);
        end
    end
    
    
    % Only use the signals withing the timeWindowForSNRanalysis
    meanSignalModulation = meanSignalModulation(timeBinsForSNR);
    noiseTraces = noiseTraces(:,timeBinsForSNR);
    
    % Concatenate all instances along a long vector
    signal = repmat(meanSignalModulation, [1 noisyInstancesNum]);
    noise  = reshape(noiseTraces', [1 numel(noiseTraces)]); 

    % Compute the SNR on these long vectors
    theSNR = measureSNR(signal, noise);

    if (displaySNRVectors)
        displayedInstancesNum = 4;
        % Just display 10 instances of the SNR signals
        displayedTimeBins = 1:(displayedInstancesNum*numel(timeBinsForSNR));

        figure(1234); 
        if (strcmp(signalName, 'cone excitations'))
            subplot(2,2,1);
            cla;
            ylabel('pRate');
        else
            subplot(2,2,2);
            cla
            ylabel('pCurrent');
        end
        stairs(signal(displayedTimeBins), 'r-', 'LineWidth', 2); hold on;
        stairs(noise(displayedTimeBins), 'k-');
        title(sprintf('SNR = %2.2f', theSNR));
        
        % Just display the SNR mean signal  
        if (strcmp(signalName, 'cone excitations'))
            subplot(2,2,3);
            cla
            ylabel('pRate');
        else
            subplot(2,2,4);
            cla
            ylabel('pCurrent');
        end
 
        stairs(timeAxis(timeBinsForSNR)*1000, meanSignalModulation, 'r-', 'LineWidth', 2);
        drawnow
    end
    
end

function theSNR = measureSNR(signal, noise)
    if (length(signal) ~= length(noise))
        error('lengths of signal and noise must match')
    end

    signalPower = sum(abs(signal).^2)/length(signal);
    noisePower = sum(abs(noise).^2)/length(signal);
    theSNR = 10 * log10(signalPower / noisePower);
end

