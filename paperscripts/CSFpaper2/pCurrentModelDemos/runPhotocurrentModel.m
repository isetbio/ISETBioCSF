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

    useDefaultImplementation = ~true;
    

    if (useDefaultImplementation)
        switch (eccentricity)
            case 'foveal'
                os = osBioPhys('eccentricity', 0);
            case 'peripheral'   
                os = osBioPhys('eccentricity', 20);
            otherwise
                error('Unknown eccentricity: ''%s''. Choose either ''foveal'' or ''peripheral''.', eccentricity);
        end
        
        % simulation time step
        dt = stimulus.timeAxis(2)-stimulus.timeAxis(1);
    
        fprintf('Computing %d response instances (default implementation)\n', noisyInstancesNum);
        os = osSet(os, 'noise flag', 'none');
        % Compute steady-state
        state = osAdaptSteadyState(os, stimulus.pRate(1));
        s = struct('state', state, 'timeStep', dt);
        % Compute mean pCurrent response to the isomerization stimuli
        pCurrents = osAdaptTemporal(stimulus.pRate, s);
        if (noisyInstancesNum>0)
            noisyPcurrents = osAddNoise(pCurrents, 'sampTime', dt);
        end
        
        % Trim all signals by removing all time points before the warmUpTime
        [~,idx0] = min(abs(stimulus.timeAxis-stimulus.warmUpTimeSeconds));

        % Keep 20 milliseconds before stimulus onset
        negativeTime = 20/1000;
        idx = find(stimulus.timeAxis >= stimulus.warmUpTimeSeconds-negativeTime);
        keptIndices = idx(1):length(stimulus.timeAxis);
    
        modelResponse.timeAxis = stimulus.timeAxis(keptIndices) - stimulus.timeAxis(idx0);
        modelResponse.pRate = stimulus.pRate(keptIndices);
        modelResponse.opsin = 0*modelResponse.timeAxis;
        modelResponse.pde = modelResponse.opsin;
        modelResponse.ca = modelResponse.opsin;
        modelResponse.gC = modelResponse.opsin;
        modelResponse.caSlow = modelResponse.opsin;
        modelResponse.cGMP = modelResponse.opsin;
        modelResponse.membraneCurrent = pCurrents(keptIndices);
        if (noisyInstancesNum>0)
            modelResponse.noisyMembraneCurrents = noisyPcurrents(keptIndices);
        end
        % Obtain adaptation current as the current at the last point in the warmup period.
        modelResponse.membraneCurrentAdaptation = pCurrents(idx(1));
        return;
    end
    
    switch (eccentricity)
        case 'foveal'
            modelConstants = struct(...
                'g_R', 12, ...           % opsin gain  
                'rho_R', 10, ...         % opsin inactivation rate constant (1/s)
                'p_dark', 700, ...       % phosphodiesterase spontaneous activation rate constant in the dark(1/s)
                'rho_PDE', 22, ...       % phosphodiesterase activity inactivation rate constant 1/s
                'n_GC',4, ...            % cooperativity for cyclase, hill coef
                'k_GC',0.5, ...          % hill affinity for cyclase
                'g_dark',20.5, ...       % concentration of cGMP in darkness
                'k_cGMP', 0.02, ...      % constant relating cGMP to current
                'n_cGMP', 3, ...         % cooperativity for cGMP->current
                'ca_dark',1, ...         % calcium concentration in darkness
                'rho_Ca', 5, ...         % rate constant for calcium removal in 1/s
                'betaSlow', 0.4 ...     % rate constant for slow calcium modulation of channels 
            );
        
        case 'peripheral'
            modelConstants = struct(...
                'g_R', 10, ...            % opsin gain  
                'rho_R', 22, ...         % opsin inactivation rate constant (1/s)
                'p_dark', 2000, ...      % phosphodiesterase spontaneous activation rate constant in the dark(1/s)
                'rho_PDE', 22, ...       % phosphodiesterase activity decay rate 1/s
                'n_GC',4, ...            % cooperativity for cyclase, hill coef
                'k_GC',0.5, ...          % hill affinity for cyclase
                'g_dark',20.5, ...       % concentration of cGMP in darkness
                'k_cGMP', 0.02, ...       % constant relating cGMP to current
                'n_cGMP', 3, ...          % cooperativity for cGMP->current
                'ca_dark',1, ...          % calcium concentration in darkness
                'rho_Ca', 9, ...          % rate constant for calcium removal in 1/s
                'betaSlow', 0.4 ...    % rate constant for slow calcium modulation of channels    
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
        pde(ii+1) = computePDEactivation(pde(ii), opsin(ii+1), modelConstants, dt);

        % Compute Ca and slow Ca concentrations
        ca(ii+1) = computeCalciumConcentration(ca(ii), Imembrane(ii), modelConstants, dt);
        caSlow(ii+1) = computeCaSlow(caSlow(ii), ca(ii+1), modelConstants, dt);
        
        % Instantaneous function of the Ca concentration
        gC(ii+1) = computeGuanlylateCyclaseActivation(ca(ii+1), modelConstants);

        % Compute cGMP concentration
        cGMP(ii+1) = computeCyclicGMPconcentration(cGMP(ii), gC(ii+1), pde(ii+1), dt);

        % compute membrane current (photocurrent) which is an instanerous function of cGMP and caSlow
        Imembrane(ii+1) = computeMembraneCurrent(cGMP(ii+1), caSlow(ii+1), modelConstants);
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
    
    % Add outer segment noise instances
    if (noisyInstancesNum > 0)
        fprintf('Computing %d response instances (components implementation - not default)\n', noisyInstancesNum);
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
    % dOpsin/dt = gR * pRate(t) - rhoR * opsin(t)
    gain = modelConstants.g_R;
    decayRateConstant = modelConstants.rho_R;
    opsin = opsin + dt * (gain * pRate - decayRateConstant * opsin);
end

function pde = computePDEactivation(pde, opsin, modelConstants, dt)
    % Implementing the following equation
    % dPDE/dt = opsin(t) + eta -phi*PDE(t)
    decayRateConstant = modelConstants.rho_PDE;
    pde = pde + dt * (opsin + modelConstants.p_dark - decayRateConstant * pde);
end

function gC = computeGuanlylateCyclaseActivation(ca, modelConstants)
    % gC is an instantaneous function of the Ca concentration
    GCmax = modelConstants.p_dark / modelConstants.rho_PDE * modelConstants.g_dark * (1 + (modelConstants.ca_dark / modelConstants.k_GC) ^ modelConstants.n_GC);
    gC = GCmax ./ (1 + (ca / modelConstants.k_GC) .^ modelConstants.n_GC);
end

function Imembrane = computeMembraneCurrent(cGMP, caSlow, modelConstants)
    Imembrane = modelConstants.k_cGMP * (cGMP^modelConstants.n_cGMP / (1+caSlow/modelConstants.ca_dark)); 
end

function ca = computeCalciumConcentration(ca, Imembrane, modelConstants, dt)
    % Ca concentration depends on calcium entry and calcium extrusion   
    q = 2 * modelConstants.rho_Ca * modelConstants.ca_dark / (modelConstants.k_cGMP * modelConstants.g_dark ^ modelConstants.n_cGMP);
    ca = ca + dt * (q * Imembrane - modelConstants.rho_Ca * ca);
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


function [adaptedData, model] = osAdaptTemporal(pRate, obj)
% Time varying current response from photon rate and initial state
%
% Syntax:
%   adaptedData = osAdaptTemporal(pRate, obj)
%
% Description:
%    This function is only called internally from @osBioPhys/osCompute.m
%
%    Time varying current response from photon rate and initial state.
%
%    In this case, the physiological differential equations for cones are
%    implemented. The differential equations are:
%
%       1) d opsin(t) / dt = -sigma * opsin(t) + R*(t)
%       2) d PDE(t) / dt = opsin(t) - phi * PDE(t) + eta
%       3) d cGMP(t) / dt = S(t) - PDE(t) * cGMP(t)
%       4) d Ca(t) / dt = q * I(t) - beta * Ca(t)
%       5) d Ca_slow(t) / dt = - beta_slow * (Ca_slow(t) - Ca(t))
%       6) S(t) = smax / (1 + (Ca(t) / kGc)^n)
%       7) I(t) = k * cGMP(t) ^ h / (1 + Ca_slow / Ca_dark)
%
%    This model gives a cone-by-cone adaptation and produces a time-series
%    structure in adaptedDat that is stored into the current field of the
%    cone mosaic object in @osBioPhys/osCompute.m.
%
%    Examples are contained in the code. To access, type 'edit
%    osAdaptTemporal.m' into the Command Window.
%
% Inputs:
%    pRate      - Photon absorption rate,
%                 coneMosaic.absorptions/coneMosaic.integrationTime.
%    obj        - osBioPhys object containing many initial parameters
%
% Outputs:
%   adaptedData - adapted photocurrent data (pA) for coneMosaic.current
%   obj         - osBioPhys object containing many final parameters
%
% Optional key/value pairs:
%    None.
%
% References:
%   http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%   https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
%
% Notes:
%    * [Note: JNM - Example doesn't work!]
%
% See Also:
%    osAdaptSteadyState, osAdaptTemporal
%

% History:
%    xx/xx/14  HJ   ISETBIO Team, 2014
%    08/xx/16  JRG  ISETBIO Team, updated 8/2016
%    02/14/18  jnm  Formatting
%    04/07/18  dhb  Skip broken example.

% Examples:
%{
    % ETTBSkip.  To work, this example will need some inputs defined before
    % the funtion is called.
    %
    % From @osBioPhys/osCompute.m, line 64:
    [current, model.state] = osAdaptTemporal(pRate, model.state);
%}

%%  Check inputs
if ~exist('pRate', 'var') || isempty(pRate)
    error('Photon absorption rate required.');
end

dt = obj.timeStep;
model = obj.state;

%% Simulate differential equations

adaptedData = zeros(1, size(pRate, 2) + 1);
adaptedData(:, 1) = model.bgCur;

q = 2 * model.beta * model.cdark / (model.k * model.gdark ^ model.h);
smax = model.eta / model.phi * model.gdark * ...
    (1 + (model.cdark / model.kGc) ^ model.n);


for ii = 1 : size(pRate, 2)
    model.opsin = model.opsin + dt * (model.OpsinGain * pRate(:, ii) ...
        - model.sigma * model.opsin);
    model.PDE = model.PDE + dt * (model.opsin + model.eta - model.phi * ...
        model.PDE);
    model.Ca = model.Ca + dt * (q * model.k * model.cGMP .^ model.h ./ ...
        (1 + model.Ca_slow / model.cdark) - model.beta * model.Ca);
    model.Ca_slow = model.Ca_slow - dt * model.betaSlow * ...
        (model.Ca_slow - model.Ca);
    model.st = smax ./ (1 + (model.Ca / model.kGc) .^ model.n);
    model.cGMP = model.cGMP  + dt * (model.st - model.PDE .* model.cGMP);

    adaptedData(:, ii) = - model.k * model.cGMP .^ model.h ./ ...
        (1 + model.Ca_slow / model.cdark);
end

adaptedData(:, size(pRate, 2) + 1) = adaptedData(:, size(pRate, 2));
adaptedData = adaptedData(:, 2:end);
% model = obj.model;
end