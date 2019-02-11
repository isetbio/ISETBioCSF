function demonstratePhotocurrentModel

% On the disk membrane absorption of a photon by an opsin pigment molecule 
% leads to instantaneous (10-13 sec) photoisomerization which transforms the
% opsin modelcule to its activated state (R*).
               
pulseDurationSeconds = 1/1000;          % 1 millisecond
stimPhotonsDeliveredDuringPulse = 1;    % 

pRatesBackground = [0 300 1000 3000 10000]; % % R*/sec
legends = {};

for backgroundIndex = 1:numel(pRatesBackground) 
    pRateBackground = pRatesBackground(backgroundIndex);
    [timeAxis, pRateStimulus, opsin, pde, ca, caSlow, cGMP, ImembraneNull] = ...
        computeAll(pRateBackground, 0, pulseDurationSeconds);

    [timeAxis, pRateStimulus, opsin, pde, ca, caSlow, cGMP, ImembraneStim] = ...
        computeAll(pRateBackground, stimPhotonsDeliveredDuringPulse, pulseDurationSeconds);

    impulseResponse(backgroundIndex,:) = ImembraneStim-ImembraneNull;
    legends{numel(legends)+1} = sprintf('%2.0f R*/cone/sec', pRateBackground);
end

figure(2); clf;
subplot(6,1,1);
plotTimeSeries(timeAxis, pRateStimulus, 'k', 'photons', 'line');
set(gca, 'YLim', [0 max(pRateStimulus)*1.4]);
subplot(6,1,2);
plotTimeSeries(timeAxis, opsin, 'k', 'R*', 'line');
subplot(6,1,3);
plotTimeSeries(timeAxis, pde, 'k', 'PDE', 'line');
subplot(6,1,4);
plotTimeSeries(timeAxis, ca, 'k', 'Ca concentration', 'line');
hold on;
plotTimeSeries(timeAxis, caSlow, 'r', 'slowCa concentration', 'line');
subplot(6,1,5);
plotTimeSeries(timeAxis, cGMP, 'k', 'cGMP concentration', 'line');
subplot(6,1,6);
plotTimeSeries(timeAxis, ImembraneStim, 'k', 'photocurrent', 'line');


plotImpulseResponses(timeAxis, impulseResponse, pRatesBackground, legends);

end

function plotImpulseResponses(timeAxis, impulseResponse, pRatesBackground, legends)
    cmap = brewermap(numel(pRatesBackground) , 'spectral');

    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1]);
    subplot(2,2,1); hold on
    for backgroundIndex = 1:numel(pRatesBackground) 
        plot((timeAxis-0.1)*1000, squeeze(impulseResponse(backgroundIndex,:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 2);
    end
    for backgroundIndex = 1:numel(pRatesBackground) 
        plot((timeAxis-0.1)*1000, squeeze(impulseResponse(backgroundIndex,:)), '-', ...
            'Color', [0 0 0], 'LineWidth', 3);
        plot((timeAxis-0.1)*1000, squeeze(impulseResponse(backgroundIndex,:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 2);
    end

    hold off
    xlabel('\it time (msec)');
    set(gca, 'YLim', max(abs(impulseResponse(:)))*[-0.3 1.1], 'XLim', [0.0 400]);
    set(gca, 'FontSize', 14);
    legend(legends);
    grid on; box on;

    subplot(2,2,3); hold on
    for backgroundIndex = 1:numel(pRatesBackground) 
        scaledIR = squeeze(impulseResponse(backgroundIndex,:));
        scaledIR = scaledIR/max(abs(scaledIR));
        plot((timeAxis-0.1)*1000, scaledIR, '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 2);
    end
    for backgroundIndex = 1:numel(pRatesBackground) 
        scaledIR = squeeze(impulseResponse(backgroundIndex,:));
        scaledIR = scaledIR/max(abs(scaledIR));
        plot((timeAxis-0.1)*1000, scaledIR, '-', ...
            'Color', [0 0 0], 'LineWidth', 3);
        plot((timeAxis-0.1)*1000, scaledIR, '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 2);
    end

    hold off
    xlabel('\it time (msec)');
    legend(legends);
    set(gca, 'YLim', [-0.3 1.1], 'XLim', [0 400]);
    set(gca, 'FontSize', 14);
    grid on; box on;

    nFFT = 4096*8;
    spectra = zeros(numel(pRatesBackground), nFFT);
    for backgroundIndex = 1:numel(pRatesBackground) 
        scaledIR = squeeze(impulseResponse(backgroundIndex,:));
        scaledIR = scaledIR/max(abs(scaledIR));
        spectra(backgroundIndex,:) = fftshift(abs(fft(scaledIR, nFFT)));
    end

    deltaT = timeAxis(2)-timeAxis(1);
    maxTF = 1/(2*deltaT);
    deltaTF = maxTF/(nFFT/2);
    tfAxis = ((-nFFT/2):(nFFT/2-1))*deltaTF;
    
    subplot(2,2,[2 4]); hold on;
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


function [timeAxis, pRate, opsin, pde, ca, caSlow, cGMP, Imembrane] = ...
    computeAll(pRateBackground, stimPhotonsDeliveredDuringPulse, pulseDurationSeconds)
%
    simulationTimeStepSeconds = 0.1/1000;
    warmUpTime = 25.0;
    stimOnset = 25.1;
    responseDurationSeconds = stimOnset + 1.0;

    % Generate time axis
    timeAxis = 0:simulationTimeStepSeconds:responseDurationSeconds;
    dt = timeAxis(2)-timeAxis(1);

    % Photon rate (photons/sec) with a single photon at stimOnsetSeconds
    pRate = zeros(1,length(timeAxis)) + pRateBackground;

    stimBins = round(pulseDurationSeconds/dt);
    stimBinIndices = round(stimOnset/dt) + (1:stimBins);

    stimPRate = stimPhotonsDeliveredDuringPulse/(pulseDurationSeconds);
    pRate(stimBinIndices) = pRate(stimBinIndices) + stimPRate;

    % Opsin activation
    opsin = zeros(1,length(timeAxis));
    pde = opsin;
    cGMP = opsin;
    gC = opsin;
    ca = opsin;
    caSlow = opsin;
    Imembrane = opsin;
    
    peripheralModelConstants = struct(...
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
        'gamma', 10 ...     % opsin gain  
    );

    fovealModelConstants = struct(...
        'sigma', 10, ...        % rhodopsin acivity decay rate constant (1/s)
        'phi', 22, ...          % phosphodiesterase activity decay rate 1/sg
        'eta', 700, ...        % phosphodiesterase spontaneous activation rate constant (1/s)
        'gdark',20.5, ...       % concentration of cGMP in darkness
        'k', 0.02, ...          % constant relating cGMP to current
        'h', 3, ...             % cooperativity for cGMP->current
        'cdark',1, ...          % calcium concentration in darkness
        'beta', 5, ...          % rate constant for calcium removal in 1/s
        'betaSlow', 0.4, ...    % rate constant for slow calcium modulation of channels
        'n',4, ...              % cooperativity for cyclase, hill coef
        'kGc',0.5, ...          % hill affinity for cyclase
        'gamma', 12 ...     % opsin gain  
    );

    modelConstants = fovealModelConstants;
    
    
    for ii = 1 : length(pRate)-1
        opsin(ii+1) = computeOpsinActivation(opsin(ii), pRate(ii), modelConstants, dt);
        pde(ii+1)   = computePDEactivation(pde(ii), opsin(ii), modelConstants, dt);

        % instantaneous function of the Ca concentration
        gC(ii)      = computeGuanlylateCyclaseActivation(ca(ii), modelConstants);

        ca(ii+1)    = computeCalciumConcentration(ca(ii), Imembrane(ii), modelConstants, dt);
        caSlow(ii+1) = computeCaSlow(caSlow(ii), ca(ii), modelConstants, dt);

        cGMP(ii+1)  = computeCyclicGMPconcentration(cGMP(ii), gC(ii), pde(ii), dt);

        % photocurrent is an instanerous function of cGMP and caSlow
        Imembrane(ii+1) = computeMembraneCurrent(cGMP(ii), caSlow(ii), modelConstants);
    end

    % Trim all signals by removing all time points before the warmUpTime
    idx = find(timeAxis > warmUpTime);
    keptIndices = idx(1):length(timeAxis);
    timeAxis = timeAxis(keptIndices) - timeAxis(idx(1));
    pRate = pRate(keptIndices);
    opsin = opsin(keptIndices);
    pde = pde(keptIndices);
    ca = ca(keptIndices);
    caSlow = caSlow(keptIndices);
    cGMP = cGMP(keptIndices);
    Imembrane =-Imembrane(keptIndices);
end


function plotTimeSeries(timeAxis, signal, lineColor, signalName, plotType)
if (strcmp(plotType, 'stem'))
    stem(timeAxis, signal, 'Color', lineColor, 'MarkerFaceColor', lineColor, 'MarkerSize', 6, 'LineWidth', 1.0);
else
    plot(timeAxis, signal, 'o', 'MarkerSize', 5, 'MarkerFaceColor', lineColor, 'Color', lineColor, 'LineWidth', 1.0);
end
    ylabel(sprintf('\\it %s',signalName));
    xlabel('\it time (seconds)');
    set(gca, 'XLim', [0 1.2], 'FontSize', 12);
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