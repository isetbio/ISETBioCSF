function testPhotocurrentModel
    recompute = true;
    if (recompute)
    cParams.spontaneousIsomerizationRate = 200; %  R*/c/s
    cParams.eccentricity  = 'foveal';
    cParams.noisyInstancesNum = 1000;
    cParams.useDefaultImplementation = true;
    
    % Define the stim params struct
    vParams.weberContrast = -95/100;
    vParams.adaptationPhotonRate = 12000;
    vParams.pulseDurationSeconds = 100/1000;
    vParams.photonIntegrationTime = vParams.pulseDurationSeconds;    % controls binning of photons in photon count signal
    
    contrastLevels = [0.02 0.04 0.08 0.16 0.32 0.64 0.90 -0.02 -0.04 -0.08 -0.16 -0.32 -0.64 -0.90 ];
    adaptationLevels = [1000 6000 20000]; % [600 1000 1800 3000 6000 10000 20000];
    
    for iAdaptationIndex = 1:numel(adaptationLevels)
        for iContrastIndex = 1:numel(contrastLevels)
            vParams.weberContrast = contrastLevels(iContrastIndex);
            vParams.adaptationPhotonRate = adaptationLevels(iAdaptationIndex);
            [~, theConeExcitationSNR(iAdaptationIndex, iContrastIndex), ...
                thePhotoCurrentSNR(iAdaptationIndex, iContrastIndex)] = runSimulation(vParams, cParams);
        end
    end
    
    save('results.mat', 'theConeExcitationSNR', 'thePhotoCurrentSNR', 'contrastLevels', 'adaptationLevels');
    else
        load('results.mat', 'theConeExcitationSNR', 'thePhotoCurrentSNR', 'contrastLevels', 'adaptationLevels');
    end
    
    % Plot results
    incrementsIdx = find(contrastLevels>0);
    decrementsIdx = find(contrastLevels<0);
    
    plotPolaritySNRs(contrastLevels(incrementsIdx), adaptationLevels, ...
        theConeExcitationSNR(:,incrementsIdx), thePhotoCurrentSNR(:,incrementsIdx), 'increments', 2);
    
    plotPolaritySNRs(-contrastLevels(decrementsIdx), adaptationLevels, ...
        theConeExcitationSNR(:,decrementsIdx), thePhotoCurrentSNR(:,decrementsIdx), 'decrements', 3);
    
end

function plotPolaritySNRs(contrastLevels, adaptationLevels, theConeExcitationSNR, thePhotoCurrentSNR, contrastPolarity, figNo)
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 300 950], 'Color', [1 1 1]);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', 1, ...
       'rowsNum', 3, ...
       'heightMargin',   0.1, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.11, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);
   
    subplot('Position', subplotPosVectors(1,1).v);
    plotSNR(contrastLevels, adaptationLevels, theConeExcitationSNR, true, true, true, true, sprintf('cone excitations (%s)', contrastPolarity));
   
    subplot('Position', subplotPosVectors(2,1).v);
    plotSNR(contrastLevels, adaptationLevels, thePhotoCurrentSNR, true, true, true, true, sprintf('photocurrent (%s)', contrastPolarity));

    theSNRratios = thePhotoCurrentSNR./theConeExcitationSNR;
    
    subplot('Position', subplotPosVectors(3,1).v);
    legends = {};
    cMap = brewermap(numel(contrastLevels), '*Spectral');
    
    for iContrastIndex = 1:numel(contrastLevels)
        color = squeeze(cMap(iContrastIndex,:));
        legends{iContrastIndex} = sprintf('%2.0f%%', contrastLevels(iContrastIndex)*100);
        plot(adaptationLevels, theSNRratios(:,iContrastIndex), 'o-', ...
            'Color', 0.5*color, 'MarkerFaceColor', color, ...
            'MarkerSize', 12, 'LineWidth', 1.5); hold on;
    end
    set(gca, 'XLim', [adaptationLevels(1)*0.9 adaptationLevels(end)*1.1], 'XTIck', [600 2000 6000 20000], 'YLim', [0.01 1.5], 'XScale', 'log');
    grid on; box on;
    axis 'square'
    set(gca, 'FontSize', 14);
    xlabel('\it adaptation level (R*/c/s)');
    ylabel(sprintf('\\it SNR ratio'));
    legend(legends, 'Location', 'NorthEast');
    title(sprintf('photocurrents : cone excitations\n(%s)',contrastPolarity));
end


function plotSNR(contrastLevels, adaptationLevels, theSNR, showXLabel, showXTickLabels, showYLabel, showYTickLabels, signalName)
   cMap = brewermap(numel(adaptationLevels), '*RdYlBu');
   legends = {};
   for iAdaptationIndex = 1:numel(adaptationLevels)
        color = squeeze(cMap(iAdaptationIndex,:));
        legends{iAdaptationIndex} = sprintf('%2.0f', adaptationLevels(iAdaptationIndex));
        plot(contrastLevels*100, theSNR(iAdaptationIndex,:), 'o-', ...
            'Color', 0.5*color, 'MarkerFaceColor', color, ...
            'MarkerSize', 12, 'LineWidth', 1.5); hold on;
    end
    set(gca, 'XLim', [0.015 1.0]*100, 'XTick', [1 3 10 30 100], 'YTick', [0.1 0.3 1 3 10 30], 'YLim', [0.1 50], 'XScale', 'log', 'YScale', 'log');
    grid on; box on;
    axis 'square'
    set(gca, 'FontSize', 14);
    
    if (showYLabel)
        ylabel('\it SNR');
    end
    if (showXLabel)
        xlabel('\it weber contrast (%)');
    end
    
    if (~showYTickLabels)
        set(gca, 'YTickLabel', {});
    end
    
    if (~showXTickLabels)
        set(gca, 'XTickLabel', {});
    end
    
    title(signalName);
    legend(legends, 'Location', 'NorthWest');
end

function [modelResponse, theConeExcitationSNR, thePhotoCurrentSNR] = runSimulation(vParams, cParams)

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
        [theConeExcitationSNR, noiseTimeOnset, coneExcitationModulationPeak, coneExcitationPhotocurrentNoiseSigma] = computeSNR(modelResponse.noisyConeExcitationTimeAxis, modelResponse.noisyConeExcitationRates);
        [thePhotoCurrentSNR, noiseTimeOnset, photocurrentModulationPeak, photocurrentNoiseSigma] = computeSNR(modelResponse.timeAxis, modelResponse.noisyMembraneCurrents);
        
        % Plot responses
        plotModelResponse(modelResponse, theConeExcitationSNR, thePhotoCurrentSNR, noiseTimeOnset, ...
            coneExcitationModulationPeak, coneExcitationPhotocurrentNoiseSigma, photocurrentModulationPeak, photocurrentNoiseSigma);
end


function [theSNR, noiseTimeOnset, modulationPeak, noiseSigma] = computeSNR(timeAxis, noisyResponses)
    % measure noise properties from the last 0.5 seconds of the signal
    noiseTimeOnset = timeAxis(end)-0.5;
    
    % Estimate noiseSigma after noiseTimeOnset
    idx = find(timeAxis>=noiseTimeOnset);
    noisyResponses2 = noisyResponses(:,idx);
    noiseSigma = std(noisyResponses2(:),  0, 1);
    
    % Estimate mean of noise 
    noiseMean = mean(noisyResponses2(:));

    % Estimate mean of signal from up to last 0.3 seconds of the signal
    idx = find(timeAxis < noiseTimeOnset);
    meanResponse = mean(noisyResponses(:,idx), 1);
    
    % Estimate peak signal modulation
    modulationPeak = max(abs(meanResponse-noiseMean));
    theSNR = modulationPeak / noiseSigma;
end



function plotModelResponse(modelResponse, theConeExcitationSNR, thePhotoCurrentSNR, noiseTimeOnset, ...
    coneExcitationModulationPeak, coneExcitationPhotocurrentNoiseSigma, photocurrentModulationPeak, photocurrentNoiseSigma)

    coneExcitationRange = [min(modelResponse.noisyConeExcitations(:))*0.9 max(modelResponse.noisyConeExcitations(:))*1.1];
    coneExcitationRateRange = [min(modelResponse.noisyConeExcitationRates(:))*0.9 max(modelResponse.noisyConeExcitationRates(:))*1.1];
    photoCurrentRange = [min(modelResponse.noisyMembraneCurrents(:)) max(modelResponse.noisyMembraneCurrents(:))];
    
    figure(1);clf;
    
    subplot(4,1,1);
    stairs(modelResponse.noisyConeExcitationTimeAxis, modelResponse.noisyConeExcitations',  'k-'); hold on;
    stairs(modelResponse.noisyConeExcitationTimeAxis, modelResponse.meanConeExcitationCountSignal, 'r-',  'LineWidth', 2.0);
    stairs(modelResponse.noisyConeExcitationTimeAxis, mean(modelResponse.noisyConeExcitations,1), 'b-',  'LineWidth', 2.0);
    plot(noiseTimeOnset*[1 1], coneExcitationRange, 'b--', 'LineWidth', 1.5);
    
    
    
    
    ylabel('\it photon absorptions (R*/c/tau)'); xlabel('\it time (sec)');
    set(gca, 'FontSize', 14);
    set(gca, 'YLim', coneExcitationRange, 'XLim', [0 modelResponse.noisyConeExcitationTimeAxis(end)]);
   
    
    subplot(4,1,2);
    stairs(modelResponse.noisyConeExcitationTimeAxis, modelResponse.meanConeExcitationRateSignal, 'r-',  'LineWidth', 2.0); hold on;
    stairs(modelResponse.noisyConeExcitationTimeAxis, mean(modelResponse.noisyConeExcitationRates,1),  'b--', 'LineWidth', 2.0); hold on;
    ylabel('\it mean photon absorption rates (R*/c/sec)');  xlabel('\it time (sec)');
    set(gca, 'XLim', [0 modelResponse.noisyConeExcitationTimeAxis(end)]);
    set(gca, 'FontSize', 14);
    title(sprintf('SNR = %2.2f', theConeExcitationSNR));
    
    subplot(4,1,3);
    stairs(modelResponse.noisyConeExcitationTimeAxis, modelResponse.noisyConeExcitationRates',  'k-'); hold on;
    stairs(modelResponse.noisyConeExcitationTimeAxis, modelResponse.meanConeExcitationRateSignal, 'r-',  'LineWidth', 2.0);
    plot(noiseTimeOnset*[1 1], coneExcitationRateRange, 'b--', 'LineWidth', 1.5);
    plot(0.5*[1 1], coneExcitationModulationPeak*[-0.5 0.5]+modelResponse.meanConeExcitationRateSignal(end), 'ms-', 'LineWidth', 2.0);
    plot(0.55*[1 1], coneExcitationPhotocurrentNoiseSigma*[-0.5 0.5]+modelResponse.meanConeExcitationRateSignal(end), 'bs-', 'LineWidth', 2.0);
    ylabel('\it photon absorption rate (R*/c/sec)');  xlabel('\it time (sec)');
    set(gca, 'YLim', coneExcitationRateRange, 'XLim', [0 modelResponse.noisyConeExcitationTimeAxis(end)]);
    set(gca, 'FontSize', 14);
    title(sprintf('SNR = %2.2f', theConeExcitationSNR));
    
    subplot(4,1,4);
    plot(modelResponse.timeAxis, modelResponse.noisyMembraneCurrents', 'k-');
    hold on;
    plot(modelResponse.timeAxis, modelResponse.membraneCurrent, 'r-',  'LineWidth', 2.0);
    plot(noiseTimeOnset*[1 1], photoCurrentRange, 'b--', 'LineWidth', 1.5);
    plot(0.5*[1 1], photocurrentModulationPeak*[-0.5 0.5] + modelResponse.membraneCurrent(end), 'ms-', 'LineWidth', 2.0);
    plot(0.55*[1 1], photocurrentNoiseSigma*[-0.5 0.5]+ modelResponse.membraneCurrent(end), 'bs-', 'LineWidth', 2.0);
    set(gca, 'YLim', photoCurrentRange, 'XLim', [0 modelResponse.noisyConeExcitationTimeAxis(end)]);
    ylabel('\it photocurrent (pA)');  xlabel('\it time (sec)');
    set(gca, 'FontSize', 14);
    title(sprintf('SNR = %2.2f', thePhotoCurrentSNR));
    
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

