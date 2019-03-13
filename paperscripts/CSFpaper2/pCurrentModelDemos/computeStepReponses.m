function [timeAxis, stepResponses, modelResponses, legends] = ...
    computeStepReponses(constantStimParams, adaptationPhotonRates, pulseWeberContrasts, pulseDurationSeconds,  spontaneousIsomerizationRate, eccentricity, noisyInstancesNum, useDefaultImplementation)
% Run the photocurrent model for a number of step stimuli
%
% Syntax:
%   [timeAxis, stepResponses, modelResponses, legends] = ...
%    computeStepReponses(constantStimParams, adaptationPhotonRates, pulseWeberContrasts, ...
%    pulseDurationSeconds, simulationTimeStepSeconds, spontaneousIsomerizationRate, eccentricity);
%
% Description:
%    Function to run the outer segment photocurrent model for a number of 
%    step stimuli, which may vary in:
%       - the adaptation rate
%       - the pulse amplitude (specified as Weber contrast)
%       - the pulse duration
%    The user can specify the simulation time step and the model
%    eccentricity.
%
% Inputs:
%    adaptationPhotonRates     - Vector of background photon rates (R*/c/s)
%    pulseWeberContrasts       - Vector of pulse Weber contrasts
%    pulseDurationSeconds      - Scalar. Pulse duration in seconds
%    spontaneousIsomerizationRate - Scalar. Spontaneous isomerization rate
%    eccentricity              - String. 'foveal' or 'peripheral'
%    noisyInstancesNum         - Number of noisy instances to compute
%    useDefaultImplementation  - true or false
%                 - true to use the os object - no internal component visualization
%                 - false to visualize the internal components
%
% Output:
%    timeAxis       - the time axis of the response
%    stepResponses  - the responses to the various step stimuli
%    modelResponses - cell array containing the various model component
%                     responses
%    legends        - cell array containing strings with the examined
%                     condition
%
% Optional key/value pairs:
%    None.

% History:
%    2/13/19  NPC   ISETBIO Team, 2019

    
    % Initialize
    stimParams = constantStimParams;
    
    nAdaptationLevels = numel(adaptationPhotonRates);
    mPulseStrengths = 2*numel(pulseWeberContrasts);
    photonsDeliveredDuringPulses = zeros(1, mPulseStrengths);
    modelResponses = cell(nAdaptationLevels, mPulseStrengths);
    legends = cell(nAdaptationLevels, mPulseStrengths);
    
    % Run model for the different adaptation levels
    for adaptationIndex = 1:nAdaptationLevels
        
        for weberIndex = 1:numel(pulseWeberContrasts)
            photonsDeliveredDuringPulse = ...
                pulseWeberContrasts(weberIndex)*adaptationPhotonRates(adaptationIndex)*pulseDurationSeconds;
            % additional photons during the increment stimulus
            photonsDeliveredDuringPulses((weberIndex-1)*2+1) = photonsDeliveredDuringPulse;    
            % additional photons during the decrement stimulus
            photonsDeliveredDuringPulses((weberIndex-1)*2+2) =-photonsDeliveredDuringPulse;
        end
    
        % Run all pulse strengths for this adaptation level
        for pulseStrengthIndex = 1:mPulseStrengths
            % Design stimulus
            stimParams.adaptationPhotonRate = adaptationPhotonRates(adaptationIndex);
            stimParams.photonsDeliveredDuringPulse = photonsDeliveredDuringPulses(pulseStrengthIndex);
            stimulus = designPhotonRateStimulus(stimParams, spontaneousIsomerizationRate);
            
            % Run model
            model = runPhotocurrentModel(stimulus, eccentricity, noisyInstancesNum, useDefaultImplementation);
            
            % Obtain step response by subtracting the adaptation membrane current
            if (adaptationIndex*pulseStrengthIndex==1)
                stepResponses = zeros(nAdaptationLevels, mPulseStrengths, length(model.membraneCurrent));
            end
            stepResponses(adaptationIndex, pulseStrengthIndex, :) = model.membraneCurrent;
        
            modelResponses{adaptationIndex, pulseStrengthIndex} = model;
            legends{adaptationIndex,pulseStrengthIndex} = ...
                sprintf('adapt:%2.0f p/c/s\npulse:%2.0f p/c/s', stimParams.adaptationPhotonRate, stimParams.adaptationPhotonRate);
            timeAxis = model.timeAxis;
        end
    end 
end