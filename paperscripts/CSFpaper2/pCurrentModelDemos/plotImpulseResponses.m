function plotImpulseResponses(timeAxis, impulseResponses, temporalFrequencyAxis, frequencySpectra, adaptationPhotonRates, legends, figNo)
% Plot the photocurrent responses for the examined impulse stimuli 
%
% Syntax:
%   plotImpulseResponses(timeAxis, impulseResponses, temporalFrequencyAxis, frequencySpectra, adaptationPhotonRates, legends, figNo)
%
% Description:
%    Plot the photocurrent impulse responses and their frequency frequencySpectra
%    for the different adaptation levels.
%    Impulse responses and frequencySpectra are plotted both unscaled and normalized
%
% Inputs:
%    timeAxis              - time axis of the response
%    impulseResponses      - matrix with impulse responses
%    temporalFrequencyAxis - temporal frequency axis
%    frequencySpectra      - matrix with frequency spectra
%    adaptationPhotonRates - vector with the examined adaptation levels
%    legends               - cell array containing strings with the examined conditions
%    figNo                 - figure number
%
% Output:
%    None.
%
% Optional key/value pairs:
%    None.

% History:
%    2/13/19  NPC   ISETBIO Team, 2019

    % Color scheme to use for the different impulse resposes
    cmap = brewermap(numel(adaptationPhotonRates) , 'spectral');

    xLims = [0 200];
    xTicks = 0:50:500;
    
    hFig = figure(figNo); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 720 720]);
    
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 3, ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.07, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.15, ...
       'topMargin',      0.1);

    subplot(2,2,1); hold on
    for backgroundIndex = 1:numel(adaptationPhotonRates) 
        plot(timeAxis*1000, squeeze(impulseResponses(backgroundIndex,:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 1.5);
    end
    for backgroundIndex = 1:numel(adaptationPhotonRates) 
        plot(timeAxis*1000, squeeze(impulseResponses(backgroundIndex,:)), '-', ...
            'Color', [0 0 0], 'LineWidth', 2);
        plot(timeAxis*1000, squeeze(impulseResponses(backgroundIndex,:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 1.5);
    end

    hold off
    xlabel('\it time (msec)');
    set(gca, 'YLim', [-0.2 0.95], 'XLim', xLims, 'XTick', xTicks);
    set(gca, 'FontSize', 14);
    legend(legends);
    grid on; box on;

    subplot(2,2,3); hold on
    for backgroundIndex = 1:numel(adaptationPhotonRates) 
        scaledIR = squeeze(impulseResponses(backgroundIndex,:));
        scaledIR = scaledIR/max(abs(scaledIR));
        plot(timeAxis*1000, scaledIR, '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 1.5);
    end
    for backgroundIndex = 1:numel(adaptationPhotonRates) 
        scaledIR = squeeze(impulseResponses(backgroundIndex,:));
        scaledIR = scaledIR/max(abs(scaledIR));
        plot(timeAxis*1000, scaledIR, '-', ...
            'Color', [0 0 0], 'LineWidth', 2);
        plot(timeAxis*1000, scaledIR, '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 1.5);
    end

    hold off
    xlabel('\it time (msec)');
    legend(legends);
    set(gca, 'YLim', [-0.22 1.02], 'XLim', xLims, 'XTick', xTicks);
    set(gca, 'FontSize', 14);
    grid on; box on;

    
    subplot(2,2,2); hold on;
    for backgroundIndex = 1:numel(adaptationPhotonRates)
        plot(temporalFrequencyAxis, squeeze(frequencySpectra(backgroundIndex,:))/max(frequencySpectra(:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 1.5);
    end
    for backgroundIndex = 1:numel(adaptationPhotonRates)
        plot(temporalFrequencyAxis, squeeze(frequencySpectra(backgroundIndex,:))/max(frequencySpectra(:)), '-', ...
            'Color', [0 0 0], 'LineWidth', 2);
        plot(temporalFrequencyAxis, squeeze(frequencySpectra(backgroundIndex,:))/max(frequencySpectra(:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 1.5);
    end
    xlabel('\it temporal frequency (Hz)');
    legend(legends, 'Location', 'SouthWest');
    tfIdx = find(temporalFrequencyAxis>=0&temporalFrequencyAxis<=100);
    set(gca, 'XLim', [temporalFrequencyAxis(tfIdx(1)) temporalFrequencyAxis(tfIdx(end))], 'XScale', 'log', 'YScale', 'log'); 
    set(gca, 'XTick', [1 3 10 30 100], 'FontSize', 14);
    grid on; box on;
    
    % Plot the normalized frequencySpectra
    frequencySpectra = bsxfun(@times, frequencySpectra, 1./frequencySpectra(:,1));
    subplot(2,2,4); hold on;
    for backgroundIndex = 1:numel(adaptationPhotonRates)
        plot(temporalFrequencyAxis, squeeze(frequencySpectra(backgroundIndex,:))/max(frequencySpectra(:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 1.5);
    end
    for backgroundIndex = 1:numel(adaptationPhotonRates)
        plot(temporalFrequencyAxis, squeeze(frequencySpectra(backgroundIndex,:))/max(frequencySpectra(:)), '-', ...
            'Color', [0 0 0], 'LineWidth', 2);
        plot(temporalFrequencyAxis, squeeze(frequencySpectra(backgroundIndex,:))/max(frequencySpectra(:)), '-', ...
            'Color', squeeze(cmap(backgroundIndex,:)), 'LineWidth', 1.5);
    end
    xlabel('\it temporal frequency (Hz)');
    legend(legends, 'Location', 'SouthWest');
    tfIdx = find(temporalFrequencyAxis>=0&temporalFrequencyAxis<=100);
    set(gca, 'XLim', [temporalFrequencyAxis(tfIdx(1)) temporalFrequencyAxis(tfIdx(end))], 'XScale', 'log', 'YScale', 'log'); 
    set(gca, 'XTick', [1 3 10 30 100], 'FontSize', 14);
    grid on; box on;    
end