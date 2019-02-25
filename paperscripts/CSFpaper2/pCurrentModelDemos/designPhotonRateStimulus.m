function stimulus = designPhotonRateStimulus(stimParams, spontaneousIsomerizationRate)
% Design a photon rate stimulus (input to runPhotocurrentModel.m)
%
% Syntax:
%   stimulus = designPhotonRateStimulus(stimParams, spontaneousIsomerizationRate);
%
% Description:
%    Function to run the outer segment photocurrent model for impulse
%    stimuli presented on different adaptation levels.
%    The user can specify the simulation time step and the model
%    eccentricity.
%
% Inputs:
%    stimParams     - struct with different parameters
%    spontaneousIsomerizationRate - the rate of spontaneous photon
%    absorptions
%
% Output:
%    stimulus       - struct with stimulus components
%
% Optional key/value pairs:
%    None.

% History:
%    2/13/19  NPC   ISETBIO Team, 2019

    % Setup warm up time
    stimulus.warmUpTimeSeconds = 25.0;
    stimulus.onsetSeconds = stimulus.warmUpTimeSeconds + stimParams.timeSampleSeconds;
    stimulus.durationSeconds = stimulus.onsetSeconds + stimParams.totalDurationSeconds;
    
    % Time axis
    stimulus.timeAxis = 0:stimParams.timeSampleSeconds:stimulus.durationSeconds;
    dt = stimParams.timeSampleSeconds;
    
    % Background
    stimulus.pRate = zeros(1,length(stimulus.timeAxis)) + stimParams.adaptationPhotonRate;
    
    % Add the spontaneous isomerization rate
    stimulus.pRate = stimulus.pRate + spontaneousIsomerizationRate;
    
    % Design stimulus
    switch (stimParams.type)
        case 'pulse'
            % Determine time bins for the pulse
            stimBins = round(stimParams.pulseDurationSeconds/dt);
            stimBinIndices = round(stimulus.onsetSeconds/dt) + (0:(stimBins-1));
            % Determine photon rate for the pulse
            pulsePhotonRate = stimParams.photonsDeliveredDuringPulse/stimParams.pulseDurationSeconds;
            % Add the pulse photon rate to the background photon rate
            stimulus.pRate(stimBinIndices) = stimulus.pRate(stimBinIndices) + pulsePhotonRate;
        otherwise
            error('Unknown stimulus type: ''%s''.', stimParams.type);
    end
end
