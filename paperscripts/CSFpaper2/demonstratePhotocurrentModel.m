function demonstratePhotocurrentModel

% On the disk membrane absorption of a photon by an opsin pigment molecule 
% leads to instantaneous (10-13 sec) photoisomerization which transforms the
% opsin modelcule to its activated state (R*).

simulationTimeStepSeconds = 0.1/1000;

stimPhotonsDeliveredDuringPulse = -100;
pulseDurationSeconds = 500/1000;
pRateBackground = 300;


stimOnset= 25.5;
warmUpTime = 25.0;
responseDurationSeconds = stimOnset + 2.0;

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

modelConstants = struct(...
    'sigma', 22, ...  % rhodopsin acivity decay rate constant (1/s)
    'gamma', 10, ...    % opsin gain
    'eta', 2000, ...     % phosphodiesterase spontaneous activation rate constant (1/s)
    'phi', 22, ...       % phosphodiesterase activity decay rate 1/sg
    'cdark',1, ...      % calcium concentration in darkness
    'gdark',20.5, ...   % concentration of cGMP in darkness
    'beta', 9, ...      % rate constant for calcium removal in 1/s
    'betaSlow', 0.4, ...  % rate constant for slow calcium modulation of channels
    'kGc',0.5, ...      % hill affinity for cyclase
    'n',4, ...           % cooperativity for cyclase, hill coef
    'k', 0.02, ...       % constant relating cGMP to current
    'h', 3 ...         % cooperativity for cGMP->current
    );

    
    
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

figure(1); clf;
subplot(6,1,1);
plotTimeSeries(timeAxis, pRate, 'k', 'photons');
set(gca, 'YLim', [0 max(pRate)*1.4]);
subplot(6,1,2);
plotTimeSeries(timeAxis, opsin, 'k', 'R*');
subplot(6,1,3);
plotTimeSeries(timeAxis, pde, 'k', 'PDE');
subplot(6,1,4);
plotTimeSeries(timeAxis, ca, 'k', 'Ca concentration');
hold on;
plotTimeSeries(timeAxis, caSlow, 'r', 'slowCa concentration');
subplot(6,1,5);
plotTimeSeries(timeAxis, cGMP, 'k', 'cGMP concentration');
subplot(6,1,6);
plotTimeSeries(timeAxis, Imembrane, 'k', 'photocurrent');
end

function plotTimeSeries(timeAxis, signal, lineColor, signalName)
    plot(timeAxis, signal, '-', 'Color', lineColor, 'LineWidth', 1.5)
    ylabel(signalName);
    xlabel('time (seconds)');
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
