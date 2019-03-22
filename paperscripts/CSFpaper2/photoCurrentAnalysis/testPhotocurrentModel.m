function testPhotocurrentModel
    
    spontaneousIsomerizationRate = 100; %  R*/c/s
    eccentricity  = 'foveal';
    noisyInstancesNum = 100;
    
    useDefaultImplementation = true;
    
    % Define the stim params struct
    weberContrast = 100/100;
    adaptationPhotonRate = 5000;
    pulseDurationSeconds = 100/1000;
    photonsDeliveredDuringPulse = weberContrast*adaptationPhotonRate*pulseDurationSeconds;
    
    photonIntegrationTime = 5/1000;    % controls binning of photons in photon count signal
    
    stimParams = struct(...
        'type', 'pulse', ...                                        % type of stimulus
        'timeSampleSeconds', 0.1/1000, ...                          % simulation time step
        'pulseDurationSeconds', pulseDurationSeconds, ...           % pulse duration in seconds
        'totalDurationSeconds', 500/1000, ...                       % total duration of the stimulus
        'adaptationPhotonRate',  adaptationPhotonRate, ...          %  R*/c/s
        'photonsDeliveredDuringPulse', photonsDeliveredDuringPulse  ...  %  R*/c during  the total pulse duration
    );

    stimulus = designPhotonRateStimulus(stimParams, spontaneousIsomerizationRate);
    modelResponse = runPhotocurrentModel(stimulus, eccentricity, noisyInstancesNum, photonIntegrationTime, useDefaultImplementation);
    plotModelResponse(modelResponse);
    
    computeResponseDetectability(modelResponse);
    
end

function computeResponseDetectability(modelResponse)
    idx = find(modelResponse.meanConeExcitationCountSignal ~= modelResponse.meanConeExcitationCountSignal(end));
    signalDistribution = modelResponse.noisyConeExcitations(:, idx);
    signalDistribution = signalDistribution(:);
  
    idx2 = setdiff(1:numel(modelResponse.noisyConeExcitationTimeAxis),  idx);
    idx2 = idx2(1:numel(idx));
    noiseDistribution = modelResponse.noisyConeExcitations(:, idx2);
    noiseDistribution = noiseDistribution(:);
    
    % Simpson & Fitter, 1973; Swets, 1986a, 1986b
    %Swets, J. A. (1986b). Indices of discrimination or diagnostic accuracy: Their ROC?s and implied models. Psychological Bulletin, 99(1), 100?117.
    %Swets, J. A. (1996). Signal Detection Theory and ROC Analysis in Psychology and Diagnostics: Collected Papers. Mahwah, NJ: Lawrence Erlbaum Associates.
    % Simpson, A. J., & Fitter, M. J. (1973). What is the best index of detectability? Psychological Bulletin, 80(6), 481?488.
    
    d_subA = (mean(signalDistribution)-mean(noiseDistribution)) / ...
                  sqrt(0.5*(var(signalDistribution)+var(noiseDistribution)));
    
    bins = linspace(min(modelResponse.noisyConeExcitations(:)), max(modelResponse.noisyConeExcitations(:)), 50);
    [countsSignal,centers] = hist(signalDistribution, bins);
    [countsNoise,centers] = hist(noiseDistribution, bins);
    figure(2)
    bar(centers,[countsSignal(:) countsNoise(:)], 1.0)
    title(sprintf('d_subA = %2.2f\n', d_subA));
end


function plotModelResponse(modelResponse)

    coneExcitationRange = [min(modelResponse.noisyConeExcitations(:))*0.9 max(modelResponse.noisyConeExcitations(:))*1.1];
    photoCurrentRange = [min(modelResponse.noisyMembraneCurrents(:)) max(modelResponse.noisyMembraneCurrents(:))];
    
    coneExcitationNoise = bsxfun(@minus, modelResponse.noisyConeExcitations, modelResponse.meanConeExcitationCountSignal);
    photoCurrentNoise = bsxfun(@minus, modelResponse.noisyMembraneCurrents, modelResponse.membraneCurrent);
    
    coneExcitationSTD = std(modelResponse.noisyConeExcitations,  0, 1);
    photoCurrentSTD = std(modelResponse.noisyMembraneCurrents, 0, 1);
    
    
    figure(1);clf;
    subplot(3,3,1);
    plot(modelResponse.timeAxis, modelResponse.pRate, 'k-');
    ylabel('photon absorption rate (R*/c/s)');
    
    subplot(3,3,2);
    stairs(modelResponse.noisyConeExcitationTimeAxis, modelResponse.noisyConeExcitations(1,:),  'k-'); hold on;
    stairs(modelResponse.noisyConeExcitationTimeAxis, mean(modelResponse.noisyConeExcitations,1), 'b-',  'LineWidth', 2.0)
    stairs(modelResponse.noisyConeExcitationTimeAxis, modelResponse.meanConeExcitationCountSignal, 'r-',  'LineWidth', 2.0)
    stairs(modelResponse.noisyConeExcitationTimeAxis, modelResponse.meanConeExcitationCountSignal + coneExcitationSTD, 'k--', 'LineWidth', 1.5);
    stairs(modelResponse.noisyConeExcitationTimeAxis, modelResponse.meanConeExcitationCountSignal - coneExcitationSTD, 'k--', 'LineWidth', 1.5);
    
    ylabel('photon absorptions (R*/c/tau)');
    set(gca, 'YLim', coneExcitationRange);
    
    subplot(3,3,3);
    stairs(modelResponse.noisyConeExcitationTimeAxis, coneExcitationNoise(1,:),  'k-');  hold  on
    stairs(modelResponse.noisyConeExcitationTimeAxis,  coneExcitationSTD,  'k--',  'LineWidth', 1.5);
    stairs(modelResponse.noisyConeExcitationTimeAxis, -coneExcitationSTD,  'k--',  'LineWidth', 1.5);
    set(gca, 'YLim', (coneExcitationRange(2)-coneExcitationRange(1))/2*[-1 1]);
    ylabel('photon absorptions (R*/c/tau)');
    
    
    subplot(3,3,5);
    plot(modelResponse.timeAxis, modelResponse.noisyMembraneCurrents(1,:), 'k-');
    hold on;
    plot(modelResponse.timeAxis, modelResponse.membraneCurrent, 'r-',  'LineWidth', 2.0);
    set(gca, 'YLim', photoCurrentRange);
    ylabel('photocurrent (pA)');
    
    subplot(3,3,6)
    plot(modelResponse.timeAxis, photoCurrentNoise(1,:), 'k-');
    set(gca, 'YLim', (photoCurrentRange(2)-photoCurrentRange(1))/2*[-1 1]);
    ylabel('photocurrent (pA)');
end