function testPhotocurrentModel

    recomputeResponses = true;
    
    % Pulse duration in seconds
    
     %vParams.pulseDurationSeconds = 50/1000;
    vParams.pulseDurationSeconds = 100/1000;
    %vParams.pulseDurationSeconds = 200/1000;
    %vParams.pulseDurationSeconds = 400/1000;
    
    % Data filename
    dataFileName = sprintf('results_%dmsec.mat', vParams.pulseDurationSeconds*1000);
        
    % Recompute responses    
    if (recomputeResponses) 
        fprintf('Will save to %s\n', dataFileName);
        
        % Constant params
        cParams.spontaneousIsomerizationRate = 200; %  R*/c/s
        cParams.eccentricity  = 'foveal';
        cParams.useDefaultImplementation = true;
        cParams.noisyInstancesNum = 2000;
        
        % Varied Params
        % Pulse contrast
        vParams.weberContrast = -95/100;
        
        % Adaptation photon rate in R*/c/s
        vParams.adaptationPhotonRate = 12000;
        
        % The photonIntegrationTime affects the SNR of the cone
        % excitations. The longer it is the higher the SNR. Setting it to 
        % the pulse duration for now
        vParams.photonIntegrationTime = vParams.pulseDurationSeconds;

        % Examined contrast levels
        contrastLevels = [-0.01 -0.02 -0.04 -0.08 -0.16 -0.32 -0.64 -1.0 0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.0];

        
        % Examined adaptation levels
        % prctile               1%      5%     10%     25%     50%     75%     90%    95%      99% 
        % adaptationLevels = [0.0317  0.1434  0.2882  0.7143  1.4192  2.1320 2.5554  2.6961  2.8136]*1e4; % for 2 deg mosaic
        adaptationLevels = [500 1000 2000 4000 8000 16000 24000]; 

        % Preallocate memory
        d = cell(numel(adaptationLevels), numel(contrastLevels));
        
        % Run all conditions
        for iAdaptationIndex = 1:numel(adaptationLevels)
            for iContrastIndex = 1:numel(contrastLevels)
                tic
                fprintf('\n Computing response %d out of %d ...', (iAdaptationIndex-1)*numel(contrastLevels) + iContrastIndex, numel(adaptationLevels)*numel(contrastLevels));
                vParams.weberContrast = contrastLevels(iContrastIndex);
                vParams.adaptationPhotonRate = adaptationLevels(iAdaptationIndex);
                d{iAdaptationIndex, iContrastIndex} = runSimulation(vParams, cParams);
                fprintf('Finished in %2.1f seconds\n', toc);
            end
        end

        save(dataFileName, 'd', 'contrastLevels', 'adaptationLevels', '-v7.3');
        fprintf('Saved data to %s\n', dataFileName);
        
    else
        % Load previously computed responses
        fprintf('Will load from %s\n', dataFileName);
        load(dataFileName, 'd', 'contrastLevels', 'adaptationLevels');       
    end
    
    for iAdaptationIndex = 1:numel(adaptationLevels)
        for iContrastIndex = 1:numel(contrastLevels)
            theConeExcitationSNR(iAdaptationIndex, iContrastIndex) = d{iAdaptationIndex, iContrastIndex}.theConeExcitationSNR;
            thePhotoCurrentSNR(iAdaptationIndex, iContrastIndex) = d{iAdaptationIndex, iContrastIndex}.thePhotoCurrentSNR;
        end
    end  
        
    if (1==2)
        % Plot results as a function of contrast
        idx = find(contrastLevels>0);
        figNo = 2;
        plotPolaritySNRs(contrastLevels(idx), adaptationLevels, ...
                theConeExcitationSNR(:,idx), thePhotoCurrentSNR(:,idx), 'increments', figNo);

        idx = find(contrastLevels<0);
        figNo = 3;
        plotPolaritySNRs(-contrastLevels(idx), adaptationLevels, ...
                theConeExcitationSNR(:,idx), thePhotoCurrentSNR(:,idx), 'decrements', figNo);
    end
    
    % Plot results as a function of adaptation photon rate
    
    % Below lims and ticks for 100 msec pulse
    SNRLims = [0.1 200];
    SNRTicks = [0.1 0.3 1 3 10 30 100 300];
    SNRratioLims = [0.01 0.8];
    SNRratioTicks = [0:0.1:1.0];
    figNo = 4;
    plotBackgroundSNR(adaptationLevels, theConeExcitationSNR, thePhotoCurrentSNR, contrastLevels,  ...
        SNRLims, SNRTicks, SNRratioLims, SNRratioTicks, vParams.pulseDurationSeconds, figNo);
    
   
    % Plot the responses
    combinePolarities = true;
    showSNR = true;
    showSNRcomponents = ~true;
    
    if (combinePolarities)
        % Plot separate responses
        for iAdaptationIndex = 1:numel(adaptationLevels)
            for iContrastIndex = 1:numel(contrastLevels)/2
                dNegativePolarity = d{iAdaptationIndex, iContrastIndex};
                dPositivePolarity = d{iAdaptationIndex, iContrastIndex+numel(contrastLevels)/2};
                figNo = 1;
                plotResponses(dNegativePolarity.modelResponse, dPositivePolarity.modelResponse,adaptationLevels(iAdaptationIndex), ...
                    contrastLevels(iContrastIndex), vParams.pulseDurationSeconds, showSNR, showSNRcomponents, figNo);
            end
        end
    else
        % Plot separate responses
        for iAdaptationIndex = 1:numel(adaptationLevels)
            for iContrastIndex = 1:numel(contrastLevels)
                dSelected = d{iAdaptationIndex, iContrastIndex};
                figNo = 1;
                plotResponses(dSelected.modelResponse, [], adaptationLevels(iAdaptationIndex), ...
                    contrastLevels(iContrastIndex), vParams.pulseDurationSeconds, showSNR, showSNRcomponents, figNo);
            end
        end
    end
    
    
    
%     plotModelResponse(dSelect.modelResponse, ...
%         dSelect.theConeExcitationSNR, ...
%         dSelect.thePhotoCurrentSNR, ...
%         dSelect.noiseEstimationLatency, ...
%         dSelect.coneExcitationModulationPeak, ...
%         dSelect.coneExcitationPhotocurrentNoiseSigma, ...
%         dSelect.photocurrentModulationPeak, ...
%         dSelect.photocurrentNoiseSigma); 
end



function plotBackgroundSNR(adaptationLevels, theConeExcitationSNR, thePhotoCurrentSNR, contrastLevels, SNRLims, SNRTicks, SNRratioLims, SNRratioTicks, pulseDurationSeconds, figNo)
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1100 460], 'Color', [1 1 1]);
    
    backgroundConeExcitationRateLims = [400 30000];
    backgroundConeExcitationRateTicks =  [300 1000 3000 10000 30000];
    backgroundConeExcitationRateTickLabels =  {'0.3k', '1k', '3k', '10k', '30k'};
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', numel(contrastLevels)/2, ...
       'rowsNum', 2, ...
       'heightMargin',   0.06, ...
       'widthMargin',    0.03, ...
       'leftMargin',     0.06, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.14, ...
       'topMargin',      0.04);

    markerSize = 10;
    for iContrastIndex = 1:(numel(contrastLevels)/2) 
        iContrastIndex2 = find(contrastLevels == -contrastLevels(iContrastIndex));
        
        color = 0.4*[1 1 1];
        subplot('Position', subplotPosVectors(1,iContrastIndex).v);
        plot(adaptationLevels, theConeExcitationSNR(:, iContrastIndex), 'ks-', ...
            'Color', 0.5*color, 'MarkerFaceColor', color, ...
            'MarkerSize', markerSize, 'LineWidth', 1.5);
        hold on;
   
        plot(adaptationLevels, thePhotoCurrentSNR(:, iContrastIndex), 'ko-', ...
            'Color', 0.5*color, 'MarkerFaceColor', color, ...
            'MarkerSize', markerSize, 'LineWidth', 1.5);
        
        color = 0.8*[1 1 1];
        plot(adaptationLevels, theConeExcitationSNR(:, iContrastIndex2), 'ks-', ...
            'Color', 0.5*color, 'MarkerFaceColor', color, ...
            'MarkerSize', markerSize, 'LineWidth', 1.5);
        
        plot(adaptationLevels, thePhotoCurrentSNR(:, iContrastIndex2), 'ko-', ...
            'Color', 0.5*color, 'MarkerFaceColor', color, ...
            'MarkerSize', markerSize, 'LineWidth', 1.5);
        
        grid on
        set(gca, 'FontSize', 14, 'XScale', 'log',  ...
            'XTick',backgroundConeExcitationRateTicks, 'XLim', backgroundConeExcitationRateLims, 'XTickLabel', {}, ...
            'YLim', SNRLims, 'YTick', SNRTicks, 'YScale', 'log');
        if (iContrastIndex == 1)
            legend({'cone exc. (decr.)', 'pCurrent (decr.)', 'cone exc. (incr.)', 'pCurrent (incr.)',}, 'Location', 'NorthWest');
            ylabel('\it SNR');
        else
            set(gca, 'YTickLabel', {});
        end
        
        title(sprintf('contrast: %2.0f%%', contrastLevels(iContrastIndex2)*100)) 
    end
    
    for iContrastIndex = 1:(numel(contrastLevels)/2) 
        legends = {};
        iContrastIndex2 = find(contrastLevels == -contrastLevels(iContrastIndex));
        
        color = 0.4*[1 1 1];
        subplot('Position', subplotPosVectors(2,iContrastIndex).v);
        plot(adaptationLevels, thePhotoCurrentSNR(:, iContrastIndex)./theConeExcitationSNR(:, iContrastIndex), 'ks-', ...
            'Color', 0.5*color, 'MarkerFaceColor', color, ...
            'MarkerSize', markerSize, 'LineWidth', 1.5);
        hold on;

        color = 0.8*[1 1 1];
        plot(adaptationLevels, thePhotoCurrentSNR(:, iContrastIndex2)./theConeExcitationSNR(:, iContrastIndex2), 'ks-', ...
            'Color', 0.5*color, 'MarkerFaceColor', color, ...
            'MarkerSize', markerSize, 'LineWidth', 1.5);

        grid on
        set(gca, 'FontSize', 14, 'XScale', 'log', ...
            'XTick',backgroundConeExcitationRateTicks, 'XLim', backgroundConeExcitationRateLims, ...
            'XTickLabel', backgroundConeExcitationRateTickLabels, ...
            'YLim', SNRratioLims, 'YTick', SNRratioTicks, 'YScale', 'linear');
        
        if (iContrastIndex == 1)
            legend({'decr.', 'incr.'}, 'Location', 'NorthEast');
            ylabel(sprintf('\\it pCurrent/cone excitations \n SNR ratio'));
            xlabel(sprintf('\\it background cone \n excitation rate (R*/c/s)'));
        else
            set(gca, 'YTickLabel', {});
        end
    end
    
    pdfFileName = sprintf('SNR_analysis_%2.1fmsec.pdf', pulseDurationSeconds*1000);
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
     
end

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
        plotModelResponse(modelResponse, theConeExcitationSNR, thePhotoCurrentSNR, noiseEstimationLatency, ...
            coneExcitationModulationPeak, coneExcitationPhotocurrentNoiseSigma, photocurrentModulationPeak, photocurrentNoiseSigma, ...
            coneExcitationPhotocurrentResponseSigma, photocurrentResponseSigma);
        
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


function [theSNR, noiseEstimationLatency, peakEstimationLatency, modulationPeak, noiseSigma, responseAtPeakSigma, responsePeak] = computeSNR(timeAxis, noisyResponses, meanResponse)
    % measure noise properties from the last 0.5 seconds of the signal
    noiseEstimationLatency = timeAxis(end)-0.5;
    
    % Estimate sigma of noise
    tBinsForNoiseEstimation = find(timeAxis>=noiseEstimationLatency);
    noisyResponses2 = noisyResponses(:,tBinsForNoiseEstimation);
    noiseSigma = std(noisyResponses2(:),  0, 1);
    
    % Estimate mean of noise 
    noiseMean = mean(noisyResponses2(:));

    % Estimate mean of signal
    tBinsForSignalEstimation = find(timeAxis < noiseEstimationLatency);
    meanResponse = meanResponse(tBinsForSignalEstimation);
    
    % Compute modulation peak
    [modulationPeak, tPeakBinIndex] = max(abs(meanResponse-noiseMean));
    responsePeak = meanResponse(tPeakBinIndex);
    
    % Estimate sigma of response at peak
    peakEstimationLatency = timeAxis(tPeakBinIndex);
    noisyResponses3 = noisyResponses(:,tPeakBinIndex);
    responseAtPeakSigma = std(noisyResponses3(:),  0, 1);
    
    % Compute the SNR
    theSNR = modulationPeak / sqrt(0.5*(noiseSigma^2 + responseAtPeakSigma^2));
    
    % Show SNR components
    showSNRcomponents = false;
    if (showSNRcomponents)
        figure(1234); clf;
        plot(timeAxis(tPeakBinIndex)*[1 1], noiseMean+[0 modulationPeak], 'rs-', 'LineWidth', 1.5); hold on
        plot(timeAxis(tPeakBinIndex)*[1 1], meanResponse(tPeakBinIndex) + responseAtPeakSigma*[-1 1], 'gs-', 'LineWidth', 1.5);
        plot(1.2*[1 1], noiseMean+noiseSigma*[-1 1], 'co-', 'LineWidth', 1.5);
        plot(timeAxis, noisyResponses(1:20,:), 'k-'); 
        plot(timeAxis(tPeakBinIndex)*[1 1], noiseMean+[0 modulationPeak], 'rs-', 'LineWidth', 1.5);
        plot(timeAxis(tPeakBinIndex)*[1 1], meanResponse(tPeakBinIndex) + responseAtPeakSigma*[-1 1], 'gs-', 'LineWidth', 1.5);
        plot(1.0*[1 1], noiseMean+noiseSigma*[-1 1], 'co-', 'LineWidth', 1.5);
        legend({'modulationPeak', 'sigmaR', 'sigmaNoise'});
        drawnow;
    end
    
end


function plotResponses(modelResponse, modelResponseOpositePolarity, adaptationLevel, contrastLevel, pulseDurationSeconds, showSNR, showSNRcomponents, figNo)
 
    if (showSNR)
        [coneExcitation.SNR, ...
         coneExcitation.noiseEstimationLatency, ...
         coneExcitation.peakEstimationLatency, ...
         coneExcitation.modulationPeak, ...
         coneExcitation.noiseSigma, ...
         coneExcitation.responseAtPeakSigma, ...
         coneExcitation.responseAtPeak] = ...
            computeSNR(modelResponse.noisyConeExcitationTimeAxis, modelResponse.noisyConeExcitations, modelResponse.meanConeExcitationCountSignal);
        
        
        [photoCurrent.SNR, ...
         photoCurrent.noiseEstimationLatency, ...
         photoCurrent.peakEstimationLatency, ...
         photoCurrent.modulationPeak, ...
         photoCurrent.noiseSigma, ...
         photoCurrent.responseAtPeakSigma, ...
         photoCurrent.responseAtPeak] = ...
            computeSNR(modelResponse.timeAxis, modelResponse.noisyMembraneCurrents, modelResponse.membraneCurrent);
        
        if (~isempty(modelResponseOpositePolarity))
            [coneExcitationOpositePolarity.SNR, ...
             coneExcitationOpositePolarity.noiseEstimationLatency, ....
             coneExcitationOpositePolarity.peakEstimationLatency, ...
             coneExcitationOpositePolarity.modulationPeak, ...
             coneExcitationOpositePolarity.noiseSigma, ...
             coneExcitationOpositePolarity.responseAtPeakSigma, ...
             coneExcitationOpositePolarity.responseAtPeak] = ...
                computeSNR(modelResponseOpositePolarity.noisyConeExcitationTimeAxis, modelResponseOpositePolarity.noisyConeExcitations, modelResponseOpositePolarity.meanConeExcitationCountSignal);
            
            [photoCurrentOpositePolarity.SNR, ...
             photoCurrentOpositePolarity.noiseEstimationLatency, ...
             photoCurrentOpositePolarity.peakEstimationLatency, ...
             photoCurrentOpositePolarity.modulationPeak, ...
             photoCurrentOpositePolarity.noiseSigma, ...
             photoCurrentOpositePolarity.responseAtPeakSigma, ...
             photoCurrentOpositePolarity.responseAtPeak] = ...
                computeSNR(modelResponseOpositePolarity.timeAxis, modelResponseOpositePolarity.noisyMembraneCurrents, modelResponseOpositePolarity.membraneCurrent);
        
        end
        
    end
    
    totalPhotonCount = pulseDurationSeconds * adaptationLevel;
    modulatedPhotonCount = totalPhotonCount * contrastLevel;
    

    switch (totalPhotonCount)
        case 100
            coneExcitationRange = [0 200];
            coneExcitationTicks = 0:25:200;
            coneExcitationTickLabels = {'0', '25', '50', '70', '100', '125', '150', '170', '200'};
            
        case 200
            coneExcitationRange = [0 400];
            coneExcitationTicks = 0:50:400;
            coneExcitationTickLabels = {'0', '50', '100', '150', '200', '250', '300', '350', '400'};
        case 400
            coneExcitationRange = [0 750];
            coneExcitationTicks = 0:100:700;
            coneExcitationTickLabels = {'0', '0.1k', '0.2k', '0.3k', '0.4k', '0.5k', '0.6k', '0.7k'};
        case 800
            coneExcitationRange = [0 1500];
            coneExcitationTicks = 0:200:1500;
            coneExcitationTickLabels = {'0', '0.2k', '0.4k', '0.6k', '0.8k', '1.0k', '1.2k', '1.4k'};
        case 1600
            coneExcitationRange = [0 3000];
            coneExcitationTicks = 0:500:3000;
            coneExcitationTickLabels = {'0', '0.5k', '1.0k', '1.5k', '2.0k', '2.5k', '3.0k', '3.5k', '4.0k', '4.5k', '5.0k', '5.5k', '6.0k'};
        case 3200
            coneExcitationRange = [0 6000];
            coneExcitationTicks = 0:500:6000;
            coneExcitationTickLabels = {'0', '0.5k', '1.0k', '1.5k', '2.0k', '2.5k', '3.0k', '3.5k', '4.0k', '4.5k', '5.0k', '5.5k', '6.0k'};
        case 6400
            coneExcitationRange = [0 11000];
            coneExcitationTicks = 0:1000:15000;
            coneExcitationTickLabels = {'0', '1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '10k', '10k', '11k'};
        otherwise
            coneExcitationRange = [0 11000];
            coneExcitationTicks = 0:1000:15000;
            coneExcitationTickLabels = {'0', '1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '10k', '10k', '11k'};
    end
    
    photoCurrentRange = [-90 0];  
    photoCurrentTicks = [-90:10:0];
    
    timeLims = [-0.1 0.8];
    timeTicks = [0:0.2:1.0];
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 350 460], 'Color', [1 1 1]);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', 1, ...
       'rowsNum', 2, ...
       'heightMargin',   0.06, ...
       'widthMargin',    0.03, ...
       'leftMargin',     0.17, ...
       'rightMargin',    0.02, ...
       'bottomMargin',   0.14, ...
       'topMargin',      0.04);
    
    extraTime = 0.5;
    dt = modelResponse.noisyConeExcitationTimeAxis(2)-modelResponse.noisyConeExcitationTimeAxis(1);
    idx = find(modelResponse.noisyConeExcitationTimeAxis>=modelResponse.noisyConeExcitationTimeAxis(end)-extraTime);
    negTime = modelResponse.noisyConeExcitationTimeAxis(idx)-modelResponse.noisyConeExcitationTimeAxis(end)-dt;
    
    % Add response component before time 0
    modelResponse.noisyConeExcitationTimeAxis = cat(2,negTime, modelResponse.noisyConeExcitationTimeAxis);
    modelResponse.noisyConeExcitations = cat(2,modelResponse.noisyConeExcitations(:,idx), modelResponse.noisyConeExcitations);
    modelResponse.meanConeExcitationCountSignal = cat(2, modelResponse.meanConeExcitationCountSignal(:,idx), modelResponse.meanConeExcitationCountSignal);
    modelResponse.noisyConeExcitationTimeAxis = modelResponse.noisyConeExcitationTimeAxis + dt;
    
    if (~isempty(modelResponseOpositePolarity))
        modelResponseOpositePolarity.noisyConeExcitationTimeAxis = cat(2,negTime, modelResponseOpositePolarity.noisyConeExcitationTimeAxis);
        modelResponseOpositePolarity.noisyConeExcitations = cat(2,modelResponseOpositePolarity.noisyConeExcitations(:,idx), modelResponseOpositePolarity.noisyConeExcitations);
        modelResponseOpositePolarity.meanConeExcitationCountSignal = cat(2, modelResponseOpositePolarity.meanConeExcitationCountSignal(:,idx), modelResponseOpositePolarity.meanConeExcitationCountSignal);
        modelResponseOpositePolarity.noisyConeExcitationTimeAxis = modelResponseOpositePolarity.noisyConeExcitationTimeAxis + dt;
    end
    
    idx = find(modelResponse.timeAxis>modelResponse.timeAxis(end)-extraTime);
    dt = modelResponse.timeAxis(2)-modelResponse.timeAxis(1);
    negTime = modelResponse.timeAxis(idx) - modelResponse.timeAxis(end)-dt;
    
    % Add response component before time0
    modelResponse.timeAxis(modelResponse.timeAxis<0) = nan;
    modelResponse.timeAxis = cat(2, negTime, modelResponse.timeAxis);
    modelResponse.noisyMembraneCurrents = cat(2, modelResponse.noisyMembraneCurrents(:,idx), modelResponse.noisyMembraneCurrents);
    modelResponse.membraneCurrent = cat(2, modelResponse.membraneCurrent(idx), modelResponse.membraneCurrent);
    if (~isempty(modelResponseOpositePolarity))
        modelResponseOpositePolarity.timeAxis(modelResponseOpositePolarity.timeAxis<0) = nan;
        modelResponseOpositePolarity.timeAxis = cat(2, negTime, modelResponseOpositePolarity.timeAxis);
        modelResponseOpositePolarity.noisyMembraneCurrents = cat(2, modelResponseOpositePolarity.noisyMembraneCurrents(:,idx), modelResponseOpositePolarity.noisyMembraneCurrents);
        modelResponseOpositePolarity.membraneCurrent = cat(2, modelResponseOpositePolarity.membraneCurrent(idx), modelResponseOpositePolarity.membraneCurrent);
    end
    
    % Display up to 500 instances
    if (isempty(modelResponseOpositePolarity))
        displayedInstancesNum = 500;
    else
        % 250 instances from one polarity, 250 from the other
        displayedInstancesNum = 250;
    end
    
    modelResponse.noisyConeExcitations = modelResponse.noisyConeExcitations(1:displayedInstancesNum,:);
    modelResponse.noisyMembraneCurrents = modelResponse.noisyMembraneCurrents(1:displayedInstancesNum,:);
    if (~isempty(modelResponseOpositePolarity))
        modelResponseOpositePolarity.noisyConeExcitations = modelResponseOpositePolarity.noisyConeExcitations(1:displayedInstancesNum,:);
        modelResponseOpositePolarity.noisyMembraneCurrents = modelResponseOpositePolarity.noisyMembraneCurrents(1:displayedInstancesNum,:);
    end

    
    subplot('Position', subplotPosVectors(1,1).v);
    plotResponseInstances = true;
    theYLabel = sprintf('\\it cone excitations (R*/c/%0.2fs)', pulseDurationSeconds);
    stairsComboPlot('g', modelResponse.noisyConeExcitationTimeAxis, modelResponse.noisyConeExcitations, modelResponse.meanConeExcitationCountSignal, ...
        coneExcitationRange, coneExcitationTicks, coneExcitationTickLabels, timeLims, timeTicks, false, true, false, theYLabel, plotResponseInstances);
    
    if (~isempty(modelResponseOpositePolarity))
        stairsComboPlot('r', modelResponseOpositePolarity.noisyConeExcitationTimeAxis, modelResponseOpositePolarity.noisyConeExcitations, ...
            modelResponseOpositePolarity.meanConeExcitationCountSignal, coneExcitationRange, coneExcitationTicks, coneExcitationTickLabels, timeLims, timeTicks, false, false, false, theYLabel, plotResponseInstances);
    end
    
    plotResponseInstances = false;
    stairsComboPlot('g', modelResponse.noisyConeExcitationTimeAxis, modelResponse.noisyConeExcitations, modelResponse.meanConeExcitationCountSignal, ...
        coneExcitationRange, coneExcitationTicks, coneExcitationTickLabels, timeLims, timeTicks, false, true, false, theYLabel, plotResponseInstances);
    if (~isempty(modelResponseOpositePolarity))
        stairsComboPlot('r', modelResponseOpositePolarity.noisyConeExcitationTimeAxis, modelResponseOpositePolarity.noisyConeExcitations, ...
            modelResponseOpositePolarity.meanConeExcitationCountSignal, coneExcitationRange, coneExcitationTicks, coneExcitationTickLabels, timeLims, timeTicks, false, false, false, theYLabel, plotResponseInstances);
    end
    
    % Show SNR components
    if (showSNR)
        if (showSNRcomponents)
            plot([0 0], modelResponse.meanConeExcitationCountSignal(end) + [0 -coneExcitation.modulationPeak], 'g-', 'LineWidth', 1.5);
            plot(coneExcitation.peakEstimationLatency+pulseDurationSeconds/2*[1 1], coneExcitation.responseAtPeak + coneExcitation.responseAtPeakSigma*[-1 1], 'gs-', 'LineWidth', 1.5); 
            plot(-0.05*[1 1], modelResponse.meanConeExcitationCountSignal(end) + coneExcitation.noiseSigma*[-1 1], 'y-', 'LineWidth', 1.5); 

            if (~isempty(modelResponseOpositePolarity))
                plot([0 0], modelResponseOpositePolarity.meanConeExcitationCountSignal(end) + [0 coneExcitationOpositePolarity.modulationPeak], 'r-', 'LineWidth', 1.5);
                plot(coneExcitationOpositePolarity.peakEstimationLatency+pulseDurationSeconds/2*[1 1], coneExcitationOpositePolarity.responseAtPeak + coneExcitationOpositePolarity.responseAtPeakSigma*[-1 1], 'rs-', 'LineWidth', 1.5); 
            end
        end
        
        if (~isempty(modelResponseOpositePolarity))
            title(sprintf('SNR(dec/inc): %2.1f/%2.1f', coneExcitation.SNR, coneExcitationOpositePolarity.SNR), 'FontWeight', 'normal');
        else
            title(sprintf('SNR: %2.1f/%2.1f', coneExcitation.SNR), 'FontWeight', 'normal');
        end
    end
    
    subplot('Position', subplotPosVectors(2,1).v);
    plotResponseInstances = true;
    linesComboPlot('g', modelResponse.timeAxis, modelResponse.noisyMembraneCurrents', modelResponse.membraneCurrent, ...
        photoCurrentRange, photoCurrentTicks, timeLims, timeTicks, false, true, true, plotResponseInstances);
    if (~isempty(modelResponseOpositePolarity))
        linesComboPlot('r', modelResponseOpositePolarity.timeAxis, modelResponseOpositePolarity.noisyMembraneCurrents', ...
            modelResponseOpositePolarity.membraneCurrent, photoCurrentRange, photoCurrentTicks, timeLims, timeTicks, false, true, true, plotResponseInstances);
    end
    
    plotResponseInstances = false;
    linesComboPlot('g', modelResponse.timeAxis, modelResponse.noisyMembraneCurrents', modelResponse.membraneCurrent, ...
        photoCurrentRange, photoCurrentTicks, timeLims, timeTicks, true, true, true, plotResponseInstances);
    if (~isempty(modelResponseOpositePolarity))
        linesComboPlot('r', modelResponseOpositePolarity.timeAxis, modelResponseOpositePolarity.noisyMembraneCurrents', ...
            modelResponseOpositePolarity.membraneCurrent, photoCurrentRange, photoCurrentTicks, timeLims, timeTicks, true, true, true, plotResponseInstances);
    end
    
    % Show SNR components
    if (showSNR)
        if (showSNRcomponents)
            plot([0 0], modelResponse.membraneCurrent(end) + [0 -photoCurrent.modulationPeak], 'g-', 'LineWidth', 1.5);
            plot(photoCurrent.peakEstimationLatency+[0 0], photoCurrent.responseAtPeak + photoCurrent.responseAtPeakSigma*[-1 1], 'gs-', 'LineWidth', 1.5); 
            plot(-0.05*[1 1], modelResponse.membraneCurrent(end) + photoCurrent.noiseSigma*[-1 1], 'y-', 'LineWidth', 1.5); 

            if (~isempty(modelResponseOpositePolarity))
                plot([0 0], modelResponseOpositePolarity.membraneCurrent(end) + [0 photoCurrentOpositePolarity.modulationPeak], 'r-', 'LineWidth', 1.5);
                plot(photoCurrentOpositePolarity.peakEstimationLatency+[0 0], photoCurrentOpositePolarity.responseAtPeak + photoCurrentOpositePolarity.responseAtPeakSigma*[-1 1], 'rs-', 'LineWidth', 1.5); 
            end
        end
    
        if (~isempty(modelResponseOpositePolarity))
            title(sprintf('SNR(dec/inc): %2.1f/%2.1f', photoCurrent.SNR, photoCurrentOpositePolarity.SNR), 'FontWeight', 'normal');
        else
            title(sprintf('SNR: %2.1f/%2.1f', photoCurrent.SNR), 'FontWeight', 'normal');
        end
    end
    
    drawnow;
    if (isempty(modelResponseOpositePolarity))
        pdfFileName = sprintf('Responses_adapt_%2.0f_contrast_%2.3f_%2.0fmsec.pdf', adaptationLevel, contrastLevel, pulseDurationSeconds*1000);
    else
        pdfFileName = sprintf('Responses_adapt_%2.0f_contrast_pm%2.3f_%2.0fmsec.pdf', adaptationLevel, abs(contrastLevel), pulseDurationSeconds*1000);
    end
    
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300)
    
end

function stairsComboPlot(color, timeAxis, responseInstances, meanResponse, responseLims, responseTicks, responseTickLabels, timeLims, timeTicks, showXLabel, showYLabel, showXTickLabel, theYLabel, plotResponseInstances)
    if (plotResponseInstances)
        [stairsTime, stairsResponse] = stairsCoords(timeAxis, responseInstances);
        hPlot = plot(stairsTime, stairsResponse,  '-', 'Color', 'k'); hold on;
        for k = 1:numel(hPlot)
            hPlot(k).Color(4) = 0.3;  % 5% transparent
        end
        hold on;
    end
    
    [stairsTime, stairsResponse] = stairsCoords(timeAxis, meanResponse);
    plot(stairsTime, stairsResponse, '-', 'Color', color,  'LineWidth', 2.0);
    if (showYLabel)
        ylabel(theYLabel);
    end
    if (showXLabel)
        xlabel('\it time (sec)');
    end
    
    grid on
    set(gca, 'FontSize', 14);
    set(gca, 'XLim', timeLims, 'XTick', timeTicks, 'YLim', responseLims, 'YTick', responseTicks, 'YTickLabel', responseTickLabels);
    
    if (~showXTickLabel)
        set(gca, 'XTickLabel', {});
    end
end

function linesComboPlot(color, timeAxis, responseInstances, meanResponse, responseRange, responseTicks, timeLims, timeTicks, showXLabel, showYLabel, showXTickLabel, plotResponseInstances)
    if (plotResponseInstances)
        hPlot = plot(timeAxis, responseInstances, '-', 'Color', 'k');
        for k = 1:numel(hPlot)
            hPlot(k).Color(4) = 0.05;  % 5% transparent
        end
        hold on;
    end
    
    plot(timeAxis, meanResponse, '-', 'Color', color,  'LineWidth', 2.0);
    grid on;
    set(gca,  'XLim', timeLims, 'XTick', timeTicks, 'YLim', responseRange, 'YTick', responseTicks);
    if (showYLabel)
        ylabel('\it photocurrent (pA)');  
    end
    if (showXLabel)
        xlabel('\it time (sec)');
    end
    if (~showXTickLabel)
        set(gca, 'XTickLabel', {});
    end
    
    set(gca, 'FontSize', 14);
end


function [xStairs, yStairs] = stairsCoords(x,y)
    xStairs(1,1) = x(1);
    yStairs(:,1) = y(:,1);
    for k = 1:numel(x)-1
        xStairs = cat(2,xStairs, x(1,k));
        yStairs = cat(2,yStairs, y(:,k));
        
        xStairs = cat(2,xStairs, x(1,k));
        yStairs = cat(2,yStairs, y(:,k+1));
    end
    xStairs = cat(2, xStairs, x(1,k));
    yStairs = cat(2, yStairs, y(:,k));
end
    
function plotModelResponse(modelResponse, theConeExcitationSNR, thePhotoCurrentSNR, noiseEstimationLatency, ...
    coneExcitationModulationPeak, coneExcitationPhotocurrentNoiseSigma, photocurrentModulationPeak, photocurrentNoiseSigma, coneExcitationPhotocurrentResponseSigma, photocurrentResponseSigma)

    coneExcitationRange = [min(modelResponse.noisyConeExcitations(:))*0.9 max(modelResponse.noisyConeExcitations(:))*1.1];
    coneExcitationRateRange = [min(modelResponse.noisyConeExcitationRates(:))*0.9 max(modelResponse.noisyConeExcitationRates(:))*1.1];
    photoCurrentRange = [min(modelResponse.noisyMembraneCurrents(:)) max(modelResponse.noisyMembraneCurrents(:))];
    
    figure(1);clf;
    
    subplot(4,1,1);
    stairs(modelResponse.noisyConeExcitationTimeAxis, modelResponse.noisyConeExcitations',  'k-'); hold on;
    stairs(modelResponse.noisyConeExcitationTimeAxis, modelResponse.meanConeExcitationCountSignal, 'r-',  'LineWidth', 2.0);
    stairs(modelResponse.noisyConeExcitationTimeAxis, mean(modelResponse.noisyConeExcitations,1), 'b-',  'LineWidth', 2.0);
    plot(noiseEstimationLatency*[1 1], coneExcitationRange, 'b--', 'LineWidth', 1.5);
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
    plot(noiseEstimationLatency*[1 1], coneExcitationRateRange, 'b--', 'LineWidth', 1.5);
    plot(0.5*[1 1], coneExcitationModulationPeak*[-0.5 0.5]+modelResponse.meanConeExcitationRateSignal(end), 'ms-', 'LineWidth', 2.0);
    plot(0.55*[1 1], coneExcitationPhotocurrentNoiseSigma*[-0.5 0.5]+modelResponse.meanConeExcitationRateSignal(end), 'bs-', 'LineWidth', 2.0);
    plot(0.6*[1 1], coneExcitationPhotocurrentResponseSigma*[-0.5 0.5]+modelResponse.meanConeExcitationRateSignal(end), 'gs-', 'LineWidth', 2.0);
    
    ylabel('\it photon absorption rate (R*/c/sec)');  xlabel('\it time (sec)');
    set(gca, 'YLim', coneExcitationRateRange, 'XLim', [0 modelResponse.noisyConeExcitationTimeAxis(end)]);
    set(gca, 'FontSize', 14);
    title(sprintf('SNR = %2.2f', theConeExcitationSNR));
    
    subplot(4,1,4);
    plot(modelResponse.timeAxis, modelResponse.noisyMembraneCurrents', 'k-');
    hold on;
    plot(modelResponse.timeAxis, modelResponse.membraneCurrent, 'r-',  'LineWidth', 2.0);
    plot(noiseEstimationLatency*[1 1], photoCurrentRange, 'b--', 'LineWidth', 1.5);
    plot(0.5*[1 1], photocurrentModulationPeak*[-0.5 0.5] + modelResponse.membraneCurrent(end), 'ms-', 'LineWidth', 2.0);
    plot(0.55*[1 1], photocurrentNoiseSigma*[-0.5 0.5]+ modelResponse.membraneCurrent(end), 'bs-', 'LineWidth', 2.0);
    plot(0.6*[1 1], photocurrentResponseSigma*[-0.5 0.5]+modelResponse.membraneCurrent(end), 'gs-', 'LineWidth', 2.0);
    set(gca, 'YLim', photoCurrentRange, 'XLim', [0 modelResponse.noisyConeExcitationTimeAxis(end)]);
    ylabel('\it photocurrent (pA)');  xlabel('\it time (sec)');
    set(gca, 'FontSize', 14);
    title(sprintf('SNR = %2.2f', thePhotoCurrentSNR));
end