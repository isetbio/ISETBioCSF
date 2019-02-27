function modelResponse = runPhotocurrentModel(stimulus, eccentricity, noisyInstancesNum)
% Run the photocurrent model for a singe stimulus
%
% Syntax:
%   modelResponse = runPhotocurrentModel(stimulus, eccentricity);
%
% Description:
%    Run the outer segment photocurrent model for a single stimulus and for
%    a specified eccentricity.
%
% Inputs:
%    stimulus     - Struct with stimulus (see designPhotonRateStimulus.m)
%    eccentricity - String. 'foveal' or 'peripheral'. Determines model
%                   constants to be used
%    noisyInstancesNum - Number of noisy instances to compute
%
% Output:
%    modelResponse - struct containing the various model component responses
%
% Optional key/value pairs:
%    None.

% History:
%    2/13/19  NPC   ISETBIO Team, 2019

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

    
    % Initialize
    opsin     = zeros(1,length(stimulus.timeAxis));
    pde       = zeros(1,length(stimulus.timeAxis));
    cGMP      = zeros(1,length(stimulus.timeAxis));
    gC        = zeros(1,length(stimulus.timeAxis));
    ca        = zeros(1,length(stimulus.timeAxis));
    caSlow    = zeros(1,length(stimulus.timeAxis));
    Imembrane = zeros(1,length(stimulus.timeAxis));
    
    % simulation time step
    dt = stimulus.timeAxis(2)-stimulus.timeAxis(1);
    
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

    % Current is negative
    Imembrane = -Imembrane;
    
    % Trim all signals by removing all time points before the warmUpTime
    [~,idx0] = min(abs(stimulus.timeAxis-stimulus.warmUpTimeSeconds));
    
    % Keep 20 milliseconds before stimulus onset
    negativeTime = 20/1000;
    idx = find(stimulus.timeAxis >= stimulus.warmUpTimeSeconds-negativeTime);
    keptIndices = idx(1):length(stimulus.timeAxis);
    
    modelResponse.timeAxis = stimulus.timeAxis(keptIndices) - stimulus.timeAxis(idx0);
    modelResponse.pRate = stimulus.pRate(keptIndices);
    modelResponse.opsin = opsin(keptIndices);
    modelResponse.pde = pde(keptIndices);
    modelResponse.ca = ca(keptIndices);
    modelResponse.gC = gC(keptIndices);
    modelResponse.caSlow = caSlow(keptIndices);
    modelResponse.cGMP = cGMP(keptIndices);
    modelResponse.membraneCurrent = Imembrane(keptIndices);
    
    % Add outer segment noise instances (32x32 cones)
    if (noisyInstancesNum > 0)
        fprintf('Computing %d response instances\n', noisyInstancesNum);
        ImembraneManyInstances = repmat(reshape(modelResponse.membraneCurrent, [1 1 numel(modelResponse.membraneCurrent)]), [1 noisyInstancesNum]);
        modelResponse.noisyMembraneCurrents = squeeze(osAddNoise(ImembraneManyInstances, 'sampTime', dt));
    else
        modelResponse.noisyMembraneCurrents = [];
    end
    
    % Obtain adaptation current as the current at the last point in the warmup period.
    modelResponse.membraneCurrentAdaptation = Imembrane(idx(1));
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
    
    cGMP = cGMP  + dt * (gC - PDE * cGMP);
end