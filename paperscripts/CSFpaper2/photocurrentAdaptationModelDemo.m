function photocurrentAdaptationModelDemo
    %old()
    
    flashDurationSeconds = 100/1000;
    flashIsomerizationsPerSec = [100 300 1000 3000 10000];
    backgroundIsomerizationsPerSec = [0 100 500 1000 2000 5000 10000 20000 50000];
    [timeAxisSeconds, isomerizationStimuli, pCurrents, noisyPcurrents] = ...
        computePcurrent(flashDurationSeconds, backgroundIsomerizationsPerSec, flashIsomerizationsPerSec);
    plotResponses(timeAxisSeconds, isomerizationStimuli, pCurrents, noisyPcurrents, backgroundIsomerizationsPerSec, flashIsomerizationsPerSec);
    
    sensitivities = computeflashSensitivities(pCurrents, backgroundIsomerizationsPerSec, flashIsomerizationsPerSec);
    plotSensitivities(sensitivities, backgroundIsomerizationsPerSec, flashIsomerizationsPerSec);
    
end

function plotSensitivities(sensitivities, backgroundIsomerizationsRates, flashIsomerizationsRates)
    nBackgrounds = numel(backgroundIsomerizationsRates);
    nflashes = numel(flashIsomerizationsRates);
    darkBackgroundIndex = 1;
    
    figure(3); clf;
    legends = {};
    for flashIndex = 1:nflashes
        legends{numel(legends)+1} = sprintf('flash: %d', flashIsomerizationsRates(flashIndex));
    end
    
    cmap = colormap(colorcube);
    for flashIndex = 1:nflashes
        darkSensitivity = sensitivities(darkBackgroundIndex,flashIndex);
        plot(backgroundIsomerizationsRates(2:end), sensitivities(2:end,flashIndex)/darkSensitivity, ...
            'rs-', 'LineWidth', 1.5, 'Color', squeeze(cmap(flashIndex,:)));
        hold on
        legend(legends)
        set(gca, 'XScale', 'log', 'Yscale', 'log', 'YLim', [0.01 1.0]);
        xlabel('background intensity (R*/sec)');
        ylabel('flash sensitivity (S/Sd)');
    end
    
end

function sensitivities = computeflashSensitivities(pCurrents, backgroundIsomerizationsRates, flashIsomerizationsRates)
    nBackgrounds = numel(backgroundIsomerizationsRates);
    nflashes = numel(flashIsomerizationsRates);
    sensitivities = zeros(nBackgrounds, nflashes);
    
    for backgroundIndex = 1:nBackgrounds
        for flashIndex = 1:nflashes
            response = squeeze(pCurrents(backgroundIndex,flashIndex,:));
            peakFlashResponseModulation = max(response) - response(1);
            flashStrength = flashIsomerizationsRates(flashIndex);
            sensitivities(backgroundIndex,  flashIndex) = peakFlashResponseModulation / flashStrength;
        end
    end
    
end


function plotResponses(timeAxisSeconds, isomerizationStimuli, pCurrents, noisyPcurrents, backgroundIsomerizationsRates, flashIsomerizationsRates)
    
    nBackgrounds = numel(backgroundIsomerizationsRates);
    nflashes = numel(flashIsomerizationsRates);
    
    hFig = figure(2); clf;
    subplotPos = NicePlot.getSubPlotPosVectors(...
       'rowsNum', nBackgrounds, ...
       'colsNum', nflashes, ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.07, ...
       'bottomMargin',   0.07, ...
       'topMargin',      0.02);
   
    for backgroundIndex = 1:nBackgrounds
        for flashIndex = 1:nflashes
            subplot('Position', subplotPos(backgroundIndex, flashIndex).v)
            plot(timeAxisSeconds, squeeze(noisyPcurrents(backgroundIndex,flashIndex,:)), 'k-', ...
                'Color', [0.4 0.4 0.4]); hold on;
            plot(timeAxisSeconds, squeeze(pCurrents(backgroundIndex,flashIndex,:)), 'c-', 'LineWidth', 3.0);
            plot(timeAxisSeconds, squeeze(pCurrents(backgroundIndex,flashIndex,:)), 'b-', 'LineWidth', 1.5);
            set(gca, 'YLim', [-95 -40], 'XLim', [-0.25 1.0]);
            title(sprintf('Bkgnd: %d, flash: %d', backgroundIsomerizationsRates(backgroundIndex), flashIsomerizationsRates(flashIndex)));
        end % flashIndex
    end % backgroundIndex
end

function [timeAxisSeconds, isomerizationStimuli, pCurrents, noisyPcurrents] = ...
       computePcurrent(flashDurationSeconds, backgroundIsomerizationsPerSec, flashIsomerizationsPerSec)
   
    warmUpTimeSeconds = 20.0;
    tSeconds = warmUpTimeSeconds + 1.5;
    tTimeStep = 0.1/1000.0;
    tBinsNum = ceil(tSeconds/tTimeStep);
    timeAxisSeconds = (0:(tBinsNum-1))*tTimeStep;
    
    flashDurationSteps = round(flashDurationSeconds/tTimeStep);
    flashDurationOnsetSeconds = warmUpTimeSeconds + 0.5;
    flashDurationOnsetSteps = round(flashDurationOnsetSeconds/tTimeStep);
    
    % Design the examined stimulus grid (backgrounds x flashes)
    nBackgrounds = numel(backgroundIsomerizationsPerSec);
    nflashes = numel(flashIsomerizationsPerSec);
    isomerizationStimuli = zeros(nBackgrounds, nflashes, tBinsNum);
    flashIndices = flashDurationOnsetSteps+(1:flashDurationSteps);
    for backgroundIndex = 1:nBackgrounds
        isomerizationStimuli(backgroundIndex, :,:) = ...
            repmat(backgroundIsomerizationsPerSec(backgroundIndex), [1 nflashes tBinsNum]);
        for flashIndex = 1:nflashes
            isomerizationStimuli(backgroundIndex,flashIndex,flashIndices) = ...
                isomerizationStimuli(backgroundIndex,flashIndex,flashIndices) + ...
                repmat(flashIsomerizationsPerSec(flashIndex), [1 1 numel(flashIndices)]);
        end
    end % backgroundIndex
    
    backgroundRates = repmat(backgroundIsomerizationsPerSec(:), [1 nflashes]);
    
    % Set up biophysics os
    os = osCreate('biophys');
    os = osSet(os, 'noise flag', 'none');
    
    % Compute steady-state
    state = osAdaptSteadyState(os, backgroundRates);
    s = struct('state', state, 'timeStep', tTimeStep);
    
    % Compute pCurrent response to isomerization stimuli
    pCurrents = osAdaptTemporal(isomerizationStimuli, s);
    
    % Add noise
    noisyPcurrents = osAddNoise(pCurrents);
    
    % Only include part of response sightly before the flashDurationOnsetSeconds
    idx = find(timeAxisSeconds>=flashDurationOnsetSeconds-0.25);
    timeAxisSeconds = timeAxisSeconds(idx)-flashDurationOnsetSeconds;
    isomerizationStimuli = isomerizationStimuli(:,:,idx);
    pCurrents = pCurrents(:,:,idx);
    noisyPcurrents = noisyPcurrents(:,:,idx);
    
end


function old
    % Cone photocurrent adaptation (desensitization)
    % From Schanpf 1990: Bright light desensitized transduction with a
    % delay.
    % In the steady state, a background of intensity I 
    % lowered the sensitivity to a weak incremental test flash by a factor 
    % 1/(1+I/Io), where Io was about 26,000 R*/sec.
    
    plotSchnapfPcurrent = false;
    if (plotSchnapfPcurrent)
        % Kinetics of photocurrent for dim flashes
        scalingFactor = 20;
        tauRiseSeconds = 25/1000;
        tauDampSeconds = 110/1000;
        tauOscillationSeconds = 220/1000;
        phaseOscillationDegs = -31;
        coeffs = [scalingFactor, tauRiseSeconds, tauDampSeconds, tauOscillationSeconds, phaseOscillationDegs];
        timeAxisSeconds = (0:1:500)/1000;
        dimFlashPcurrentResponse = coneEmpiricalDimFlash(coeffs, timeAxisSeconds);
        figure(1); clf;
        plot(timeAxisSeconds, dimFlashPcurrentResponse, 'r-');
    end
    
    
    warmUpTimeSeconds = 20.0;
    tSeconds = warmUpTimeSeconds + 2;
    tTimeStep = 0.1/1000.0;
    tBinsNum = ceil(tSeconds/tTimeStep);
    timeAxisSeconds = (0:(tBinsNum-1))*tTimeStep;
    
    flashDurationSeconds = 100/1000;
    flashDurationSteps = round(flashDurationSeconds/tTimeStep);
    flashDurationOnsetSeconds = warmUpTimeSeconds + 0.5;
    flashDurationOnsetSteps = round(flashDurationOnsetSeconds/tTimeStep);
    
    % Simulate a 2x2 mosaic with different isomerizations
    backgroundIsomerizationsPerSec = [1.78 1.78*2^4; 1.78*2^8 1.78*2^12];
    pulseIsomerizationsPerSec = 1000;
    
    isomerizationRate = repmat(backgroundIsomerizationsPerSec, [1 1 tBinsNum]);
    isomerizationRate(:,:,flashDurationOnsetSteps+(1:flashDurationSteps)) = ...
        isomerizationRate(:,:,flashDurationOnsetSteps+(1:flashDurationSteps)) + repmat(pulseIsomerizationsPerSec, [1 1 flashDurationSteps]);
    
    % Set up biophysics os
    os = osCreate('biophys');
    os = osSet(os, 'noise flag', 'none');
    
    % Compute steady-state
    state = osAdaptSteadyState(os, backgroundIsomerizationsPerSec);
    s = struct('state', state, 'timeStep', tTimeStep);
    
    % Compute pCurrent response to isomerization stimuli
    pCurrent = osAdaptTemporal(isomerizationRate, s);
    
    % Add noise
    noisyPcurrent = osAddNoise(pCurrent);
    
    
    % Only include part of response sightly before the flashDurationOnsetSeconds
    idx = find(timeAxisSeconds>=flashDurationOnsetSeconds-0.5);
    timeAxisSeconds = timeAxisSeconds(idx)-flashDurationOnsetSeconds;
    isomerizationRate = isomerizationRate(:,:,idx);
    pCurrent = pCurrent(:,:,idx);
    noisyPcurrent = noisyPcurrent(:,:,idx);
    
    subplotPos = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 4, ...
       'colsNum', 3, ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.07, ...
       'bottomMargin',   0.07, ...
       'topMargin',      0.02);

    figure(2); clf;
    subplot('Position', subplotPos(1,1).v);
    plot(timeAxisSeconds, squeeze(isomerizationRate(1,1,:)), 'k-');
    set(gca, 'YLim', backgroundIsomerizationsPerSec(1,1)+[-100 pulseIsomerizationsPerSec+100]);
    subplot('Position', subplotPos(2,1).v);
    plot(timeAxisSeconds, squeeze(isomerizationRate(1,2,:)), 'r-');
    set(gca, 'YLim', backgroundIsomerizationsPerSec(1,2)+[-100 pulseIsomerizationsPerSec+100]);
    subplot('Position', subplotPos(3,1).v);
    plot(timeAxisSeconds, squeeze(isomerizationRate(2,1,:)), 'b-');
    set(gca, 'YLim', backgroundIsomerizationsPerSec(2,1)+[-100 pulseIsomerizationsPerSec+100]);
    subplot('Position', subplotPos(4,1).v);
    plot(timeAxisSeconds, squeeze(isomerizationRate(2,2,:)), 'm-');
    set(gca, 'YLim', backgroundIsomerizationsPerSec(2,2)+[-100 pulseIsomerizationsPerSec+100]);
    
    subplot('Position', subplotPos(1,2).v);
    plot(timeAxisSeconds, squeeze(pCurrent(1,1,:)), 'k-');
    set(gca, 'YLim', [-100 -50]);
    subplot('Position', subplotPos(2,2).v);
    plot(timeAxisSeconds, squeeze(pCurrent(1,2,:)), 'r-');
    set(gca, 'YLim', [-100 -50]);
    subplot('Position', subplotPos(3,2).v);
    plot(timeAxisSeconds, squeeze(pCurrent(2,1,:)), 'b-');
    set(gca, 'YLim', [-100 -50]);
    subplot('Position', subplotPos(4,2).v);
    plot(timeAxisSeconds, squeeze(pCurrent(2,2,:)), 'm-');
    set(gca, 'YLim', [-100 -50]);
    
    subplot('Position', subplotPos(1,3).v);
    plot(timeAxisSeconds, squeeze(noisyPcurrent(1,1,:)), 'k-');
    set(gca, 'YLim', [-100 -50]);
    subplot('Position', subplotPos(2,3).v);
    plot(timeAxisSeconds, squeeze(noisyPcurrent(1,2,:)), 'r-');
    set(gca, 'YLim', [-100 -50]);
    subplot('Position', subplotPos(3,3).v);
    plot(timeAxisSeconds, squeeze(noisyPcurrent(2,1,:)), 'b-');
    set(gca, 'YLim', [-100 -50]);
    subplot('Position', subplotPos(4,3).v);
    plot(timeAxisSeconds, squeeze(noisyPcurrent(2,2,:)), 'm-');
    set(gca, 'YLim', [-100 -50]);
    
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
adaptedData = zeros([size(model.opsin) size(pRate, 3) + 1]);
adaptedData(:, :, 1) = model.bgCur;

q = 2 * model.beta * model.cdark / (model.k * model.gdark ^ model.h);
smax = model.eta / model.phi * model.gdark * ...
    (1 + (model.cdark / model.kGc) ^ model.n);

for ii = 1 : size(pRate, 3)
    model.opsin = model.opsin + dt * (model.OpsinGain * pRate(:, :, ii) ...
        - model.sigma * model.opsin);
    model.PDE = model.PDE + dt * (model.opsin + model.eta - model.phi * ...
        model.PDE);
    model.Ca = model.Ca + dt * (q * model.k * model.cGMP .^ model.h ./ ...
        (1 + model.Ca_slow / model.cdark) - model.beta * model.Ca);
    model.Ca_slow = model.Ca_slow - dt * model.betaSlow * ...
        (model.Ca_slow - model.Ca);
    model.st = smax ./ (1 + (model.Ca / model.kGc) .^ model.n);
    model.cGMP = model.cGMP  + dt * (model.st - model.PDE .* model.cGMP);

    adaptedData(:, :, ii) = - model.k * model.cGMP .^ model.h ./ ...
        (1 + model.Ca_slow / model.cdark);
end

adaptedData(:, :, size(pRate, 3) + 1) = adaptedData(:, :, size(pRate, 3));
adaptedData = adaptedData(:, :, 2:end);
% model = obj.model;
end

