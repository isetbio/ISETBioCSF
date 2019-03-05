function [timeAxis, impulseResponses, temporalFrequencyAxis, impulseResponseSpectra, ...
    modelResponses, legends] = computeImpulseReponses(impulseDurationSeconds, photonCountDuringImpulse, adaptationPhotonRates, ...
    simulationTimeStepSeconds, eccentricity)
% Run the photocurrent model for a number of impulse stimuli
%
% Syntax:
%   [timeAxis, impulseResponses, temporalFrequencyAxis, impulseResponseSpectra, ...
%    modelResponses, legends] = computeStepReponses(adaptationPhotonRates, ...
%    pulseWeberContrasts, pulseDurationSeconds, simulationTimeStepSeconds, eccentricity);
%
% Description:
%    Function to run the outer segment photocurrent model for impulse
%    stimuli presented on different adaptation levels.
%    The user can specify the simulation time step and the model
%    eccentricity.
%
% Inputs:
%    adaptationPhotonRates     - Vector of background photon rates (R*/c/s)
%    simulationTimeStepSeconds - Scalar. Simulation time step in seconds
%    eccentricity              - String. 'foveal' or 'peripheral'
%
% Output:
%    timeAxis              - the time axis of the response
%    impulseResponses      - the responses to the various impulse stimuli
%    temporalFrequencyAxis - the frequency axis for the frequency spectra
%    impulseResponseSpectra- the frequency spectra of the impulse responses
%    modelResponses        - cell array containing the various model component
%                            responses
%    legends               - cell array containing strings with the examined
%                            condition
%
% Optional key/value pairs:
%    None.

% History:
%    2/13/19  NPC   ISETBIO Team, 2019

    % Define the stim params struct
    stimParams = struct(...
        'type', 'pulse', ...                            % type of stimulus
        'adaptationPhotonRate', [], ...               % background pRate
        'pulseDurationSeconds', impulseDurationSeconds, ...             % pulse duration in seconds
        'photonsDeliveredDuringPulse', photonCountDuringImpulse, ...           % how many photons during the pulse duration
        'totalDurationSeconds', 1000/1000, ...                  % total duration of the stimulus
        'timeSampleSeconds', simulationTimeStepSeconds ...
    );

      
    % Initialize
    modelResponses = cell(1,numel(adaptationPhotonRates));
    legends = cell(1, numel(adaptationPhotonRates));
    noisyInstancesNum = 0;
    
    % Run model for the different adaptation levels
    for adaptationIndex = 1:numel(adaptationPhotonRates) 
        % Design stimulus
        stimParams.adaptationPhotonRate = adaptationPhotonRates(adaptationIndex);
        stimulus = designPhotonRateStimulus(stimParams, 0);
        
        % Run model
        model = runPhotocurrentModel(stimulus, eccentricity, noisyInstancesNum);

        % Obtain impulse response by subtracting the adaptation membrane current
        ir = model.membraneCurrent-model.membraneCurrentAdaptation;
        
        % Compute spectrum of impulse response
        deltaT = model.timeAxis(2)-model.timeAxis(1);
        maxTF = 1/(2*deltaT);
        deltaTF = 0.435;
        nFFT = round(2.0*maxTF/deltaTF);
        if (mod(nFFT,2) == 1)
            nFFT = nFFT+1;
        end
        
        % Zero pad 
        signal = zeros(1, nFFT);
        margin = nFFT - length(ir);
        signal(round(margin/2)+(1:length(ir))) = ir;
        % Compute the FFT
        irSpectrum = fftshift(abs(fft(signal, nFFT)));
        
        if (adaptationIndex == 1)
            % compute temporal frequency axis support
            temporalFrequencyAxis = ((-nFFT/2):(nFFT/2-1))*deltaTF;
            idx = find((temporalFrequencyAxis>=0)&(temporalFrequencyAxis<300));
            temporalFrequencyAxis = temporalFrequencyAxis(idx);
            impulseResponses = zeros(numel(adaptationPhotonRates), length(ir));
            impulseResponseSpectra = zeros(numel(adaptationPhotonRates), length(temporalFrequencyAxis));
        end
        
        modelResponses{adaptationIndex} = model;
        impulseResponses(adaptationIndex,:) = ir;
        impulseResponseSpectra(adaptationIndex,:) = irSpectrum(idx);
        legends{adaptationIndex} = sprintf('bkgnd: %2.0f photons/cone/s', stimParams.adaptationPhotonRate);
        timeAxis = model.timeAxis;
    end
end
