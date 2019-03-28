function dStruct = runSimulation(vParams, cParams)

        % Assemble full set of stimParams
        photonsDeliveredDuringPulse = vParams.weberContrast * vParams.adaptationPhotonRate * vParams.pulseDurationSeconds;
        
        stimParams = struct(...
            'type', 'pulse', ...                                            % type of stimulus
            'timeSampleSeconds', 0.1/1000, ...                              % simulation time step
            'pulseDurationSeconds', vParams.pulseDurationSeconds, ...       % pulse duration in seconds
            'totalDurationSeconds', vParams.pulseDurationSeconds + 1.3, ... % total duration of the stimulus
            'adaptationPhotonRate',  vParams.adaptationPhotonRate, ...      %  R*/c/s
            'photonsDeliveredDuringPulse', photonsDeliveredDuringPulse  ... %  R*/c during  the total pulse duration
        );

        % Design stimulus
        stimulus = designStimulus(stimParams, cParams.spontaneousIsomerizationRate);
        
        % Run model
        modelResponse = photocurrentModel(stimulus, cParams.eccentricity, cParams.noisyInstancesNum, vParams.photonIntegrationTime, cParams.useDefaultImplementation);
        
        % Compute SNRs
        [theConeExcitationSNR, noiseEstimationLatency, coneExcitationPeakEstimationLatency, coneExcitationModulationPeak, coneExcitationPhotocurrentNoiseSigma, coneExcitationPhotocurrentResponseSigma, coneExcitationPhotocurrentResponsePeak] = ...
            computeSNR(modelResponse.noisyConeExcitationTimeAxis, modelResponse.noisyConeExcitations, modelResponse.meanConeExcitationCountSignal);
        [thePhotoCurrentSNR, noiseEstimationLatency, photocurrentPeakEstimationLatency, photocurrentModulationPeak, photocurrentNoiseSigma, photocurrentResponseSigma, photocurrentResponsePeak] = ...
            computeSNR(modelResponse.timeAxis, modelResponse.noisyMembraneCurrents, modelResponse.membraneCurrent);
        
        % Plot responses
        plotModel = false;
        if (plotModel)
            plotModelResponse(modelResponse, theConeExcitationSNR, thePhotoCurrentSNR, noiseEstimationLatency, ...
                coneExcitationModulationPeak, coneExcitationPhotocurrentNoiseSigma, photocurrentModulationPeak, ...
                photocurrentNoiseSigma, coneExcitationPhotocurrentResponseSigma, photocurrentResponseSigma);
        end
        
        % To save space, save single precision data
        fnames = fieldnames(modelResponse);
        for fk = 1:numel(fnames)
            eval(sprintf('modelResponse.%s = single(modelResponse.%s);', fnames{fk}, fnames{fk}));
        end
        
        % Return results struct
        dStruct = struct(....
            'modelResponse', modelResponse, ...
        	'theConeExcitationSNR', theConeExcitationSNR, ...
        	'thePhotoCurrentSNR', thePhotoCurrentSNR, ...
        	'coneExcitationModulationPeak', coneExcitationModulationPeak, ...
            'coneExcitationPhotocurrentNoiseSigma', coneExcitationPhotocurrentNoiseSigma, ...
        	'photocurrentModulationPeak', photocurrentModulationPeak, ...
        	'photocurrentNoiseSigma', photocurrentNoiseSigma, ...
            'noiseEstimationLatency', noiseEstimationLatency);
end


