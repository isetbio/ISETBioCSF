function demonstratePhotocurrentModel
    figNo = 0;
    eccentricity = 'foveal';
    simulationTimeStepSeconds = 0.1/1000;
    
    %figNo = demoImpulseResponses(eccentricity, simulationTimeStepSeconds, figNo);
    
    figNo = demoStepResponses(eccentricity, simulationTimeStepSeconds, figNo);
end

function figNo = demoImpulseResponses(eccentricity, simulationTimeStepSeconds, figNo)
    % Examined adaptation levels (photons/cone/sec)
    adaptationPhotonRates = [0 300 1000 3000 10000];

    % Compute the impulse response at different adaptation levels
    [timeAxis, impulseResponses, temporalFrequencyAxis, impulseResponseSpectra, modelResponses, legends] = ...
        computeImpulseReponses(adaptationPhotonRates, simulationTimeStepSeconds, eccentricity);
    
    % Plot the impulse response at different adaptation levels
    figNo = figNo + 1;
    plotImpulseResponses(timeAxis, impulseResponses, temporalFrequencyAxis, impulseResponseSpectra, adaptationPhotonRates, legends, figNo);
    
    % Plot the model response to the impulse stimuli
    figNo = figNo + 1;
    plotModelResponses(modelResponses, legends,figNo);
end


function figNo = demoStepResponses(eccentricity, simulationTimeStepSeconds, figNo)
    % Compute the step responses at different adaptation levels
    adaptationPhotonRates = [100 300 1000 3000 10000];
    pulseDurationSeconds = 100/1000;
    pulseWeberContrasts = [0.125 0.25 0.5 0.75 1.0];
    
    [timeAxis, stepResponses, modelResponses, legends] = computeStepReponses(adaptationPhotonRates, pulseWeberContrasts, pulseDurationSeconds, simulationTimeStepSeconds, eccentricity);
    
    % Plot the different step responses
    figNo = figNo + 1;
    plotStepResponses(timeAxis, stepResponses, adaptationPhotonRates, pulseWeberContrasts, legends, figNo);
end

function [timeAxis, stepResponses, modelResponses, legends] = computeStepReponses(adaptationPhotonRates, pulseWeberContrasts, pulseDurationSeconds, simulationTimeStepSeconds, eccentricity)
    % Define stim params struct
    stimParams = struct(...
        'type', 'pulse', ...                                    % type of stimulus
        'adaptationPhotonRate', [], ...                         % background pRate
        'pulseDurationSeconds', pulseDurationSeconds, ...       % pulse duration in seconds
        'photonsDeliveredDuringPulse', [], ...                  % how many photons during the pulse duration
        'totalDurationSeconds', 0.6, ...                        % total duration of the stimulus
        'timeSampleSeconds', simulationTimeStepSeconds ...
    );

    
    % Initialize
    nAdaptationLevels = numel(adaptationPhotonRates);
    mPulseStrengths = 2*numel(pulseWeberContrasts);  
    modelResponses = cell(nAdaptationLevels, mPulseStrengths);
    legends = cell(nAdaptationLevels, mPulseStrengths);
    
    % Run model for the different adaptation levels
    for adaptationIndex = 1:nAdaptationLevels
        
        photonsDeliveredDuringPulse = [];
        for weberIndex = 1:numel(pulseWeberContrasts)
            photonsDeliveredDuringPulse(numel(photonsDeliveredDuringPulse)+1) = ...
                pulseWeberContrasts(weberIndex)*adaptationPhotonRates(adaptationIndex)*pulseDurationSeconds;
            photonsDeliveredDuringPulse(numel(photonsDeliveredDuringPulse)+1) = ...
                -pulseWeberContrasts(weberIndex)*adaptationPhotonRates(adaptationIndex)*pulseDurationSeconds;
        end
    
        for pulseStrengthIndex = 1:mPulseStrengths
            % Design stimulus
            stimParams.adaptationPhotonRate = adaptationPhotonRates(adaptationIndex);
            stimParams.photonsDeliveredDuringPulse = photonsDeliveredDuringPulse(pulseStrengthIndex);
            stimulus = designPhotonRateStimulus(stimParams);
            
            % Run model
            model = runPhotocurrentModel(stimulus, eccentricity);
            
            % Obtain step response by subtracting the adaptation membrane current
            if (adaptationIndex*pulseStrengthIndex==1)
                stepResponses = zeros(nAdaptationLevels, mPulseStrengths, length(model.membraneCurrent));
            end
            stepResponses(adaptationIndex, pulseStrengthIndex, :) = model.membraneCurrent;
        
            modelResponses{adaptationIndex} = model;
            legends{adaptationIndex,pulseStrengthIndex} = sprintf('adapt:%2.0f p/c/s\npulse:%2.0f p/c/s', stimParams.adaptationPhotonRate, stimParams.adaptationPhotonRate);
            timeAxis = model.timeAxis;
        end
    end

    
end

function [timeAxis, impulseResponses, temporalFrequencyAxis, impulseResponseSpectra, modelResponses, legends] = computeImpulseReponses(adaptationPhotonRates, simulationTimeStepSeconds, eccentricity)
    % Define stim params struct
    stimParams = struct(...
        'type', 'pulse', ...                            % type of stimulus
        'adaptationPhotonRate', 2000, ...               % background pRate
        'pulseDurationSeconds', 1/1000, ...             % pulse duration in seconds
        'photonsDeliveredDuringPulse', 1, ...           % how many photons during the pulse duration
        'totalDurationSeconds', 0.6, ...                  % total duration of the stimulus
        'timeSampleSeconds', simulationTimeStepSeconds ...
    );
      
    % Initialize
    legends = {};
    modelResponses = cell(1,numel(adaptationPhotonRates));
    
    % Run model for the different adaptation levels
    for adaptationIndex = 1:numel(adaptationPhotonRates) 
        % Design stimulus
        stimParams.adaptationPhotonRate = adaptationPhotonRates(adaptationIndex);
        stimulus = designPhotonRateStimulus(stimParams);
        
        % Run model
        model = runPhotocurrentModel(stimulus, eccentricity);

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
        legends{numel(legends)+1} = sprintf('%2.0f photons/cone/sec', stimParams.adaptationPhotonRate);
        timeAxis = model.timeAxis;
    end
end


function plotModelResponses(modelResponses, legends, figNo)
    nModels = length(modelResponses);
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1]);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', nModels, ...
       'rowsNum', 6, ...
       'heightMargin',   0.03, ...
       'widthMargin',    0.06, ...
       'leftMargin',     0.07, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);


    for iModel = 1:nModels
        model = modelResponses{iModel};

        if (iModel == 1)
            labelYaxis = true;
        else
            labelYaxis = false;
        end
        
        subplot('Position', subplotPosVectors(1,iModel).v);
        plotTimeSeries(model.timeAxis, model.pRate, 'r', sprintf('absorbed photon rate\n(photons/cone/sec)'), 'line', false, labelYaxis);
        yLim = model.pRate(1)+[0 1100];
        yTicks = 0:200:100000;
        set(gca, 'YLim', yLim, 'YTick', yTicks, 'YTickLabel', sprintf('%4.0f\n',yTicks));
        grid on; box on;
        title(legends{iModel}, 'FontSize', 12);
        
        subplot('Position', subplotPosVectors(2,iModel).v);
        plotTimeSeries(model.timeAxis, model.opsin, 'r', 'opsin activation', 'line', false, labelYaxis);
        yLim = model.opsin(1)+[0 13];
        yTicks = 0:2:20000;
        set(gca, 'YLim', yLim, 'YTick', yTicks, 'YTickLabel', sprintf('%4.0f\n',yTicks));
        grid on; box on;
        
        subplot('Position', subplotPosVectors(3,iModel).v);
        plotTimeSeries(model.timeAxis, model.pde, 'r', 'PDE', 'line', false, labelYaxis);
        yLim = model.pde(1)+[-0.05 0.4];
        yTicks = 30:0.1:600;
        set(gca, 'YLim', yLim, 'YTick', yTicks);
        grid on; box on;
        
        subplot('Position', subplotPosVectors(4,iModel).v);
        plotTimeSeries(model.timeAxis, model.ca, 'r', '', 'line', false, labelYaxis);
        hold on;
        plotTimeSeries(model.timeAxis, model.caSlow, 'b', 'Ca & slowCa', 'line', false, labelYaxis);
        yLim = model.ca(1)+[-0.003 0.001];
        yTicks = 0.004:0.001:1;
        set(gca, 'YLim', yLim, 'YTick', yTicks);
        grid on; box on;
        
        subplot('Position', subplotPosVectors(5,iModel).v);
        plotTimeSeries(model.timeAxis, model.cGMP, 'r', 'cGMP', 'line', false, labelYaxis);
        yTicks = 13:0.02:21;
        yLim = model.cGMP(1)+[-0.08 0.02];
        set(gca, 'YLim', yLim, 'YTick', yTicks, 'YTickLabel', sprintf('%2.2f\n', yTicks));
        grid on; box on;
        
        subplot('Position', subplotPosVectors(6,iModel).v);
        plotTimeSeries(model.timeAxis, model.membraneCurrent, 'r', sprintf('photocurrent\n(pAmps)'), 'line', true, labelYaxis);
        yTicks = -90:0.2:0;
        yLim = model.membraneCurrent(1)+[-0.3 1];
        set(gca, 'YLim', yLim, 'YTick', yTicks, 'YTickLabel', sprintf('%2.1f\n', yTicks));
        grid on; box on;
    end
end

function plotTimeSeries(timeAxis, signal, lineColor, signalName, plotType, labelXaxis, labelYaxis)
    if (strcmp(plotType, 'stem'))
        stem(timeAxis*1000, signal, 'Color', lineColor, 'MarkerFaceColor', lineColor, 'MarkerSize', 6, 'LineWidth', 1.0);
    elseif (strcmp(plotType, 'line'))
        plot(timeAxis*1000, signal, '-', 'Color', lineColor, 'LineWidth', 1.5);
    elseif (strcmp(plotType, 'dashed line'))
        plot(timeAxis*1000, signal, '--', 'Color', lineColor, 'LineWidth', 1.5);
    elseif (strcmp(plotType, 'dotted line'))
        plot(timeAxis*1000, signal, ':', 'Color', lineColor, 'LineWidth', 1.5);
    else
        error('Uknown plotType ''%s'' in plotTimeSeries', plotType);
    end
    if (labelYaxis)
        ylabel(sprintf('\\it %s',signalName));
    end
    
    xTicks = [0 100 200 300 400 500 600 700 800 900 1000];
    if (labelXaxis)
        xlabel('\it time (msec)');
        xTickLabels = {'0', '', '200', '', '400', '', '600', '', '800', '', '1000'};
    else
        xTickLabels = {};
    end
    set(gca, 'XLim', [timeAxis(1) timeAxis(end)]*1000, 'XTick', xTicks, 'XTickLabel', xTickLabels, 'FontSize', 14);
end

function plotStepResponses(timeAxis, stepResponses, adaptationPhotonRates, pulseWeberContrasts, legends, figNo)

    nAdaptationLevels = size(stepResponses,1);
    mPulseStrengths = size(stepResponses,2);
    
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1]);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', nAdaptationLevels, ...
       'rowsNum', mPulseStrengths/2, ...
       'heightMargin',   0.03, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.07, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);
   
   for pulseStrengthIndex = 1:2:mPulseStrengths
        labelXaxis = false;
        if (pulseStrengthIndex == mPulseStrengths-1)
           labelXaxis = true;
        end
           
        for adaptationIndex = 1:nAdaptationLevels
           labelYaxis = false;
           if (adaptationIndex == 1)
               labelYaxis = true;
           end
           weberContrastIndex = floor(pulseStrengthIndex/2)+1;
           subplot('Position', subplotPosVectors(weberContrastIndex, adaptationIndex).v);
           % step increment response
           responseIncrement = squeeze(stepResponses(adaptationIndex, pulseStrengthIndex,:));
           % step decrement response
           responseDecrement = squeeze(stepResponses(adaptationIndex, pulseStrengthIndex+1,:));
           plotTimeSeries(timeAxis, responseIncrement, 'r', '', 'line', labelXaxis, labelYaxis);
           hold on;
           plotTimeSeries(timeAxis, responseDecrement, 'b', '', 'line', labelXaxis, labelYaxis);
           plotTimeSeries(timeAxis, responseDecrement(1)-(responseDecrement-responseDecrement(1)), 'b', 'current (pAmps)', 'dotted line', labelXaxis, labelYaxis);
           text(150, -5, sprintf('%2.0f photons/cone/sec', pulseWeberContrasts(weberContrastIndex)*adaptationPhotonRates(adaptationIndex)), 'FontSize', 12);
           yLim = [-90 0];
           set(gca, 'YLim',  yLim, 'YTick', [-90:10:0]);
           if (adaptationIndex > 1)
               set(gca, 'YTickLabel', {});
           end
           grid on; box on 
           if (pulseStrengthIndex == 1)
               title(sprintf('%2.0f photons/cone/sec', adaptationPhotonRates(adaptationIndex)), 'FontSize', 12);
           end
       end
   end
   
    
end

function plotImpulseResponses(timeAxis, impulseResponse, tfAxis, spectra, pRatesBackground, legends, figNo)
    cmap = brewermap(numel(pRatesBackground) , 'spectral');

    xLims = [0 500];
    xTicks = 0:100:500;
    
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1]);
    subplot(2,2,1); hold on
    for backgroundIndex = 1:numel(pRatesBackground) 
        plot(timeAxis*1000, squeeze(impulseResponse(backgroundIndex,:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 2);
    end
    for backgroundIndex = 1:numel(pRatesBackground) 
        plot(timeAxis*1000, squeeze(impulseResponse(backgroundIndex,:)), '-', ...
            'Color', [0 0 0], 'LineWidth', 3);
        plot(timeAxis*1000, squeeze(impulseResponse(backgroundIndex,:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 2);
    end

    hold off
    xlabel('\it time (msec)');
    set(gca, 'YLim', max(abs(impulseResponse(:)))*[-0.3 1.1], 'XLim', xLims, 'XTick', xTicks);
    set(gca, 'FontSize', 14);
    legend(legends);
    grid on; box on;

    subplot(2,2,3); hold on
    for backgroundIndex = 1:numel(pRatesBackground) 
        scaledIR = squeeze(impulseResponse(backgroundIndex,:));
        scaledIR = scaledIR/max(abs(scaledIR));
        plot(timeAxis*1000, scaledIR, '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 2);
    end
    for backgroundIndex = 1:numel(pRatesBackground) 
        scaledIR = squeeze(impulseResponse(backgroundIndex,:));
        scaledIR = scaledIR/max(abs(scaledIR));
        plot(timeAxis*1000, scaledIR, '-', ...
            'Color', [0 0 0], 'LineWidth', 3);
        plot(timeAxis*1000, scaledIR, '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 2);
    end

    hold off
    xlabel('\it time (msec)');
    legend(legends);
    set(gca, 'YLim', [-0.3 1.1], 'XLim', xLims, 'XTick', xTicks);
    set(gca, 'FontSize', 14);
    grid on; box on;

    
    subplot(2,2,2); hold on;
    for backgroundIndex = 1:numel(pRatesBackground)
        plot(tfAxis, squeeze(spectra(backgroundIndex,:))/max(spectra(:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 2);
    end
    for backgroundIndex = 1:numel(pRatesBackground)
        plot(tfAxis, squeeze(spectra(backgroundIndex,:))/max(spectra(:)), '-', ...
            'Color', [0 0 0], 'LineWidth', 3);
        plot(tfAxis, squeeze(spectra(backgroundIndex,:))/max(spectra(:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 2);
    end
    xlabel('\it temporal frequency (Hz)');
    legend(legends, 'Location', 'SouthWest');
    tfIdx = find(tfAxis>=0&tfAxis<=100);
    set(gca, 'XLim', [tfAxis(tfIdx(1)) tfAxis(tfIdx(end))], 'XScale', 'log', 'YScale', 'log'); 
    set(gca, 'XTick', [1 3 10 30 100], 'FontSize', 14);
    grid on; box on;
    
    % Plot the normalized spectra
    spectra = bsxfun(@times, spectra, 1./spectra(:,1));
    subplot(2,2,4); hold on;
    for backgroundIndex = 1:numel(pRatesBackground)
        plot(tfAxis, squeeze(spectra(backgroundIndex,:))/max(spectra(:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 2);
    end
    for backgroundIndex = 1:numel(pRatesBackground)
        plot(tfAxis, squeeze(spectra(backgroundIndex,:))/max(spectra(:)), '-', ...
            'Color', [0 0 0], 'LineWidth', 3);
        plot(tfAxis, squeeze(spectra(backgroundIndex,:))/max(spectra(:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 2);
    end
    xlabel('\it temporal frequency (Hz)');
    legend(legends, 'Location', 'SouthWest');
    tfIdx = find(tfAxis>=0&tfAxis<=100);
    set(gca, 'XLim', [tfAxis(tfIdx(1)) tfAxis(tfIdx(end))], 'XScale', 'log', 'YScale', 'log'); 
    set(gca, 'XTick', [1 3 10 30 100], 'FontSize', 14);
    grid on; box on;    
end


function stimulus = designPhotonRateStimulus(stimParams)
    % Setup warm up time
    stimulus.warmUpTimeSeconds = 25.0;
    stimulus.onsetSeconds = stimulus.warmUpTimeSeconds + stimParams.timeSampleSeconds;
    stimulus.durationSeconds = stimulus.onsetSeconds + stimParams.totalDurationSeconds;
    stimulus.timeAxis = 0:stimParams.timeSampleSeconds:stimulus.durationSeconds;
    dt = stimParams.timeSampleSeconds;
    
    % Background
    stimulus.pRate = zeros(1,length(stimulus.timeAxis)) + stimParams.adaptationPhotonRate;
    
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

function model = runPhotocurrentModel(stimulus, eccentricity)
%
    dt = stimulus.timeAxis(2)-stimulus.timeAxis(1);

    % Initialize
    opsin = zeros(1,length(stimulus.timeAxis));
    pde = opsin;
    cGMP = opsin;
    gC = opsin;
    ca = opsin;
    caSlow = opsin;
    Imembrane = opsin;
    
    switch (eccentricity)
        case 'foveal'
            modelConstants = struct(...
                'sigma', 10, ...        % rhodopsin acivity decay rate constant (1/s)
                'phi', 22, ...          % phosphodiesterase activity decay rate 1/sg
                'eta', 700, ...         % phosphodiesterase spontaneous activation rate constant (1/s)
                'gdark',20.5, ...       % concentration of cGMP in darkness
                'k', 0.02, ...          % constant relating cGMP to current
                'h', 3, ...             % cooperativity for cGMP->current
                'cdark',1, ...          % calcium concentration in darkness
                'beta', 5, ...          % rate constant for calcium removal in 1/s
                'betaSlow', 0.4, ...    % rate constant for slow calcium modulation of channels
                'n',4, ...              % cooperativity for cyclase, hill coef
                'kGc',0.5, ...          % hill affinity for cyclase
                'gamma', 12 ...         % opsin gain  
            );
        
        case 'peripheral'
            modelConstants = struct(...
                'sigma', 22, ...        % rhodopsin acivity decay rate constant (1/s)
                'phi', 22, ...          % phosphodiesterase activity decay rate 1/sg
                'eta', 2000, ...        % phosphodiesterase spontaneous activation rate constant (1/s)
                'gdark',20.5, ...       % concentration of cGMP in darkness
                'k', 0.02, ...          % constant relating cGMP to current
                'h', 3, ...             % cooperativity for cGMP->current
                'cdark',1, ...          % calcium concentration in darkness
                'beta', 9, ...          % rate constant for calcium removal in 1/s
                'betaSlow', 0.4, ...    % rate constant for slow calcium modulation of channels
                'n',4, ...              % cooperativity for cyclase, hill coef
                'kGc',0.5, ...          % hill affinity for cyclase
                'gamma', 10 ...         % opsin gain  
            );
        otherwise
            error('Unknown eccentricity: ''%s''. Choose either ''foveal'' or ''peripheral''.', eccentricity)
    end
    

    % run model
    for ii = 1:length(stimulus.pRate)-1
        % compute opsin and PDE activity
        opsin(ii+1) = computeOpsinActivation(opsin(ii), stimulus.pRate(ii), modelConstants, dt);
        pde(ii+1) = computePDEactivation(pde(ii), opsin(ii), modelConstants, dt);

        % Instantaneous function of the Ca concentration
        gC(ii) = computeGuanlylateCyclaseActivation(ca(ii), modelConstants);

        % Compute Ca and slow Ca concentrations
        ca(ii+1) = computeCalciumConcentration(ca(ii), Imembrane(ii), modelConstants, dt);
        caSlow(ii+1) = computeCaSlow(caSlow(ii), ca(ii), modelConstants, dt);

        % Compute cGMP concentration
        cGMP(ii+1) = computeCyclicGMPconcentration(cGMP(ii), gC(ii), pde(ii), dt);

        % compute membrane current (photocurrent) which is an instanerous function of cGMP and caSlow
        Imembrane(ii+1) = computeMembraneCurrent(cGMP(ii), caSlow(ii), modelConstants);
    end

    % Trim all signals by removing all time points before the warmUpTime
    [~,idx0] = min(abs(stimulus.timeAxis-stimulus.warmUpTimeSeconds));
    
    % Keep 20 milliseconds before stimulus onset
    negativeTime = 20/1000;
    idx = find(stimulus.timeAxis >= stimulus.warmUpTimeSeconds-negativeTime);
    keptIndices = idx(1):length(stimulus.timeAxis);
    
    model.timeAxis = stimulus.timeAxis(keptIndices) - stimulus.timeAxis(idx0);
    model.pRate = stimulus.pRate(keptIndices);
    model.opsin = opsin(keptIndices);
    model.pde = pde(keptIndices);
    model.ca = ca(keptIndices);
    model.caSlow = caSlow(keptIndices);
    model.cGMP = cGMP(keptIndices);
    model.membraneCurrent =-Imembrane(keptIndices);
    
    % Obtain adaptation current as the current at the last point in the warmup period.
    model.membraneCurrentAdaptation = -Imembrane(idx(1));
end

function opsin = computeOpsinActivation(opsin, pRate, modelConstants, dt)
    % Implementing the following equation
    % dOpsin/dt = gamma * pRate(t) - sigma * opsin(t)
    opsin = opsin + dt * (modelConstants.gamma * pRate - modelConstants.sigma * opsin);
end

function pde = computePDEactivation(pde, opsin, modelConstants, dt)
    % Implementing the following equation
    % dPDE/dt = opsin(t) + eta -phi*PDE(t)
    pde = pde + dt * (opsin + modelConstants.eta - modelConstants.phi * pde);
end

function gC = computeGuanlylateCyclaseActivation(ca, modelConstants)
    % gC is an instantaneous function of the Ca concentration
    smax = modelConstants.eta / modelConstants.phi * modelConstants.gdark * (1 + (modelConstants.cdark / modelConstants.kGc) ^ modelConstants.n);
    gC = smax ./ (1 + (ca / modelConstants.kGc) .^ modelConstants.n);
    
end

function Imembrane = computeMembraneCurrent(cGMP, caSlow, modelConstants)
    Imembrane = modelConstants.k * (cGMP^modelConstants.h / (1+caSlow/modelConstants.cdark)); 
end

function ca = computeCalciumConcentration(ca, Imembrane, modelConstants, dt)
    % Ca concentration depends on calcium entry and calcium extrusion   
    q = 2 * modelConstants.beta * modelConstants.cdark / (modelConstants.k * modelConstants.gdark ^ modelConstants.h);
    ca = ca + dt * (q * Imembrane - modelConstants.beta * ca);
end

function caSlow = computeCaSlow(caSlow, ca, modelConstants, dt)
    % Implementing the following equation
    % d caSlow(t) / dt = -betaSlow*(caSlow(t)-ca(t))
    caSlow = caSlow - dt * modelConstants.betaSlow * (caSlow - ca);
end

function cGMP = computeCyclicGMPconcentration(cGMP, gC, PDE, dt)
    % cGMP is dictated by a balance of rate of creation by guanlylate
    % cyclase, gC, and destruction by PDE
    % Implementing the following equation
    % d cGMP(t) / dt = gC(t) - PDE(t) * cGMP(t)
    
    % The ST depends on the 
    cGMP = cGMP  + dt * (gC - PDE * cGMP);
end