function visualizeStabilizedAndDynamicsSpectra(sData, sDataStabilized, figNo)

    % Limits
    sfLims = [0 30];
    tfLims = [-100 100];
    cLims = [-20 55]; %  [-60 45]; % in dB
    
    hFig = figure(figNo+100); clf;
    plotSummarySlices(sDataStabilized, sData, cLims);
    pause
    
    hFig = figure(figNo+100); clf;
    set(hFig, 'Position', [10 10 2250 625]);
    plotSlices(sDataStabilized,'STABILIZED', sfLims, tfLims, cLims,0);
    plotSlices(sData,'DYNAMIC', sfLims, tfLims, cLims,1);
    drawnow;
end

function plotSummarySlices(sDataStabilized, sDataDynamic,  cLims)

    % organization is [TF, SFy, SFx]
    % Get the slice corresponding to sfY  0 c/deg
    targetSFy = 0;
    [~,sfYIndex] = min(abs(sDataDynamic.spatialFrequencySupport)-targetSFy);
    XTspectraStabilized = squeeze(sDataStabilized.meanSpatioTemporalSpectalDensity(:,sfYIndex,:));
    XTspectraDynamic = squeeze(sDataDynamic.meanSpatioTemporalSpectalDensity(:,sfYIndex,:));
    
    
    % Average along 4 quadrants
    [averageXTspectraDynamic, sfSupport, tfSupport] = ...
        averageAcrossQuadrantsSpectum(sDataDynamic.spatialFrequencySupport, sDataDynamic.tfSupport, XTspectraDynamic);
    [averageXTspectraStabilized, ~,~] = ...
        averageAcrossQuadrantsSpectum(sDataDynamic.spatialFrequencySupport, sDataDynamic.tfSupport, XTspectraStabilized);
    
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 4, ...
       'heightMargin',  0.1, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.04, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.08, ...
       'topMargin',      0.03);

    sfLims = [1 80];
    tfLims = [tfSupport(1) 100];
    [X,Y] = meshgrid(sfSupport, tfSupport);
    
    figure(22); clf;
    
    % The spatiotemporal spectrum of the stabilized stimulus
    ax = subplot('Position', subplotPosVectors(1,1).v);
    dbSpectrum = 10*log10(averageXTspectraStabilized);
    dbSpectrum(dbSpectrum<cLims(1)) = cLims(1);
    dbSpectrum(dbSpectrum>cLims(2)) = cLims(2);
    [c, h] = contourf(X,Y,dbSpectrum, [cLims(1):cLims(end)]);
    set(h,'LineColor','none')
    xlabel('spatial frequency, X (c/deg)');
    ylabel('temporal frequency (Hz)');
    title('stabilized stimulus');
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', sfLims, 'YLim', tfLims, 'CLim', cLims, 'ZLim', cLims);
    set(gca, 'XTick', [0.1 1 3 10 30 100], 'YTick', [0.1 1 3 10 30 100], 'FontSize', 14);
    %axis 'square';
    cBar = colorbar('horizontal', 'Location', 'northoutside');
    cBar.Label.String = 'power (dB)';
    colormap(ax, hot(1024));
    
    % The spatiotemporal spectrum of the dynamic stimulus
    ax = subplot('Position', subplotPosVectors(1,2).v);
    dbSpectrum = 10*log10(averageXTspectraDynamic);
    dbSpectrum(dbSpectrum<cLims(1)) = cLims(1);
    dbSpectrum(dbSpectrum>cLims(2)) = cLims(2);
    [c, h] = contourf(X,Y,dbSpectrum, [cLims(1):cLims(end)]);
    set(h,'LineColor','none')
    xlabel('spatial frequency, X (c/deg)');
    ylabel('temporal frequency (Hz)');
    title('dynamic stimulus');
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', sfLims, 'YLim', tfLims, 'CLim', cLims, 'ZLim', cLims);
    set(gca, 'XTick', [0.1 1 3 10 30 100], 'YTick', [0.1 1 3 10 30 100], 'FontSize', 14);
    %axis 'square';
    cBar = colorbar('horizontal', 'Location', 'northoutside');
    cBar.Label.String = 'power (dB)';
    colormap(ax, hot(1024));

    
    
    
    
    sampledTFs = [0 2 5 7 10 20 50 100];
    sampledTFsNum = numel(sampledTFs);
    
    sampledSFs = [0 3 10 20 40 60];
    sampledSFsNum = numel(sampledSFs);
    
    lineColors = brewermap(max([sampledSFsNum sampledTFsNum]) , 'Spectral');
    
    legendsTF = {};
    for k = 1:sampledTFsNum
        [~, tfIndex(k)] = min(abs(tfSupport-sampledTFs(k)));
        diffEnergyAsAFunctionOfSF(k,:) =  10*log10(averageXTspectraDynamic(tfIndex(k),:)) - 10*log10(averageXTspectraStabilized(tfIndex(k),:));
        legendsTF{k} = sprintf('%2.0fHz', tfSupport(tfIndex(k)));
    end
    
    legendsSF = {};
    for k = 1:sampledSFsNum
        [~, sfIndex(k)] = min(abs(sfSupport-sampledSFs(k)));
        diffEnergyAsAFunctionOfTF(k,:) =  10*log10(averageXTspectraDynamic(:,sfIndex(k))) - 10*log10(averageXTspectraStabilized(:,sfIndex(k)));
        legendsSF{k} = sprintf('%2.0fc/deg', sfSupport(sfIndex(k)));
    end
    
    
    lowSFs = find(sfSupport<=5);
    stimSFs = find(sfSupport>10 & sfSupport < 12);
    tfs = find(tfSupport > sfSupport(1) & tfSupport < 100)
    
    noisePowerStabilized = sum(sum(averageXTspectraStabilized(tfs, lowSFs)))
    stimulusPowerStabilized = sum(sum(averageXTspectraStabilized(tfs, stimSFs)))
    snrStabilized =  stimulusPowerStabilized/noisePowerStabilized
    
    sum(averageXTspectraStabilized(:)) 
    sum(averageXTspectraDynamic(:))
    pause
    
    subplot('Position', subplotPosVectors(1,3).v); hold on
    plot(sfSupport, sum(averageXTspectraStabilized,1), 'k-');
    hold on;
    plot(sfSupport, sum(averageXTspectraDynamic,1), 'r-');
    hL = legend({'stabilized', 'dynamic'},  'Location', 'northeast');
    set(gca, 'XLim', sfLims, 'XScale',  'log', 'FontSize', 14);
    set(gca, 'XTick', [1 3 10 30 60 100]);
    xlabel('spatial frequency, X (c/deg)');
    ylabel('power (dB)');
    box on; grid on
    
    % Spectral slices at different TFs of the stabilized stimulus
    subplot('Position', subplotPosVectors(2,1).v); hold on
    for k = 1:sampledTFsNum
%        if (k == 3)
        plot(sfSupport, 10*log10(averageXTspectraStabilized(tfIndex(k),:)), 'k-', 'Color', squeeze(lineColors(k,:)), 'LineWidth', 2);
%        end
    end
    for k = 1:sampledTFsNum
%        if (k == 3)
        plot(sfSupport, 10*log10(averageXTspectraStabilized(tfIndex(k),:)), 'k-', 'Color', squeeze(lineColors(k,:))*0.7, 'LineWidth', 4);
%        end
    end
    for k = 1:sampledTFsNum
%        if (k == 3)
        plot(sfSupport, 10*log10(averageXTspectraStabilized(tfIndex(k),:)), 'k-', 'Color', squeeze(lineColors(k,:)), 'LineWidth', 2);
%        end
    end
    %axis 'square'
    title('stabilized stimulus');
    box on; grid on;
    hL = legend(legendsTF, 'NumColumns',2, 'Location', 'northeast');
    set(gca, 'XLim', sfLims, 'YLim', cLims, 'XScale', 'log', 'FontSize', 14);
    set(gca, 'XTick', [1 3 10 30 60 100]);
    xlabel('spatial frequency, X (c/deg)');
    ylabel('power (dB)');
    
    
    
    noisePowerDynamic = sum(sum(averageXTspectraDynamic(tfs, lowSFs)))
    stimulusPowerDynamic = sum(sum(averageXTspectraDynamic(tfs, stimSFs)))
    snrDynamic =  stimulusPowerDynamic/noisePowerDynamic
    
    % Spectral (sf) slices at different TFs of the dynamic stimulus
    subplot('Position', subplotPosVectors(2,2).v); hold on
    for k = 1:sampledTFsNum
%        if (k == 3)
        plot(sfSupport, 10*log10(averageXTspectraDynamic(tfIndex(k),:)), 'k-', 'Color', squeeze(lineColors(k,:)), 'LineWidth', 2);
%        end
    end
    for k = 1:sampledTFsNum
%        if (k == 3)
        plot(sfSupport, 10*log10(averageXTspectraDynamic(tfIndex(k),:)), 'k-', 'Color', squeeze(lineColors(k,:))*0.7, 'LineWidth', 4);
%        end
    end
    for k = 1:sampledTFsNum
%        if (k == 3)
        plot(sfSupport, 10*log10(averageXTspectraDynamic(tfIndex(k),:)), 'k-', 'Color', squeeze(lineColors(k,:)), 'LineWidth', 2);
%        end
    end
    %axis 'square'
    title('dynamic stimulus');
    box on; grid on;
    hL = legend(legendsTF, 'NumColumns',2, 'Location', 'northeast');
    set(gca, 'XLim', sfLims, 'YLim', cLims, 'XScale', 'log', 'FontSize', 14);
    set(gca, 'XTick', [1 3 10 30 60 100]);
    xlabel('spatial frequency, X (c/deg)');
    ylabel('power (dB)');

    
    % Spectral (sf) slices of differential power (dynamic-stabilized) stimulus
    subplot('Position', subplotPosVectors(2,3).v); hold on
    % integrate over TF
    tfBands = find(tfSupport >= sfSupport(1) & tfSupport < 100);
    totalDynamicPower = sum(averageXTspectraDynamic(tfBands,:),1);
    totalStabilizedPower = sum(averageXTspectraStabilized,1);
    plot(sfSupport, 10*log10(totalDynamicPower), 'r-', 'LineWidth', 1.5);
    hold on
    plot(sfSupport, 10*log10(totalStabilizedPower), 'k-', 'LineWidth', 1.5);
    hL = legend({'dynamic', 'stabilized'},  'Location', 'northeast');
    set(gca, 'XLim', sfLims, 'YLim', cLims, 'XScale',  'log', 'FontSize', 14);
    set(gca, 'XTick', [1 3 10 30 60 100]);
    xlabel('spatial frequency, X (c/deg)');
    ylabel('power (dB)');
    box on; grid on
    
    
    
    
    
end

function [averageSpectra, sfAxis, tfAxis] = averageAcrossQuadrantsSpectum(sfSupport, tfSupport, xtSpectra)
    targetSFx = 0;
    [~,sfIndex0] = min(abs(sfSupport)-targetSFx);
    sfIndices = (sfIndex0+1):numel(sfSupport);
    targetTF = 0;
    [~, tfIndex0] = min(abs(tfSupport)-targetTF);
    tfIndices = (tfIndex0+1):numel(tfSupport);
    
    averageSpectra = zeros(numel(tfIndices)+1, numel(sfIndices)+1);
    % The sf = 0, tf = 0 slices
    averageSpectra(1,:) = xtSpectra(tfIndex0,sfIndex0:end);
    averageSpectra(:,1) = xtSpectra(tfIndex0:end,sfIndex0);
    
    [~,idx] = max(xtSpectra(:));
    [midRow, midCol] = ind2sub(size(xtSpectra), idx)
    [tfIndex0, sfIndex0 ]
    size(tfSupport) 
    size(sfSupport)
    figure()
    imagesc(xtSpectra(midRow+(-3:3), midCol+(-3:3)))
    size(xtSpectra)
    axis 'image';
    drawnow;
    pause
    
    % average the remaining slices across 4 quadrants
    for sf = 1:numel(sfIndices)
        sfIndex = sfIndices(sf);
        sfIndex2 = sfIndex0 - (sfIndex-sfIndex0);
        for tf = 1:numel(tfIndices)
            tfIndex = tfIndices(tf);
            tfIndex2 = tfIndex0 - (tfIndex-tfIndex0);
            averageSpectra(tf+1,sf+1) = 0.25*(...
                xtSpectra(tfIndex, sfIndex) + ...
                xtSpectra(tfIndex2, sfIndex) + ...
                xtSpectra(tfIndex2, sfIndex2) + ...
                xtSpectra(tfIndex, sfIndex2));
        end
    end
    
    sum(xtSpectra(:))
    sum(averageSpectra(:))
    pause
    
    sfIndices = sfIndex0:numel(sfSupport);
    tfIndices = tfIndex0:numel(tfSupport);
    
    sfAxis = sfSupport(sfIndices);
    tfAxis = tfSupport(tfIndices);
    
    % no zero
    deltaSF = sfAxis(2)-sfAxis(1);
    deltaTF = tfAxis(2)-tfAxis(1);
    sfAxis(1) = sfAxis(2)-deltaSF/2;
    tfAxis(1) = tfAxis(2)-deltaTF/2;
    
end

function plotSlices(sData, dataLabel, sfLims, tfLims, cLims, row)

    [maxPower,idx] = max(sData.meanSpatioTemporalSpectalDensity(:));
    [maxTFindex, maxSFYindex, maxSFXindex] = ind2sub(size(sData.meanSpatioTemporalSpectalDensity),idx);
    fprintf('Max PSD is %2.2f dB at %2.1f Hz, %2.1f c/deg (x), %2.1f c/deg (y)\n', ...
        maxPower,sData.tfSupport(maxTFindex), sData.spatialFrequencySupport(maxSFXindex), sData.spatialFrequencySupport(maxSFYindex));
    fprintf('Total power (%s): %g\n', dataLabel, 10*log10(sum(sData.meanSpatioTemporalSpectalDensity(:))));
    
    % Spatiotemporal spectra at peak SFY
    XTspectra = squeeze(sData.meanSpatioTemporalSpectalDensity(:,maxSFYindex,:))';
    % Spatiotemporal spectra at peak SFX
    YTspectra = squeeze(sData.meanSpatioTemporalSpectalDensity(:,:,maxSFXindex))';
    % Spatiotemporal spectra at peak TFX
    XYspectra = squeeze(sData.meanSpatioTemporalSpectalDensity(maxTFindex,:,:));
    
    sampledSFsNum = 15;
    sampledTFsNum = 15;
    
    lineColors = brewermap(sampledSFsNum , 'Spectral');
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 6, ...
       'heightMargin',  0.10, ...
       'widthMargin',   0.03, ...
       'leftMargin',    0.02, ...
       'rightMargin',   0.03, ...
       'bottomMargin',  0.10, ...
       'topMargin',     0.03);

    subplot('Position', subplotPosVectors(1+row,1).v);
    imagesc(sData.spatialFrequencySupport, sData.spatialFrequencySupport, 10*log10(XYspectra));
    axis 'square';
    set(gca, 'XLim', max(sfLims)*[-1/2 1/2], 'YLim', sfLims , 'FontSize', 14);
    set(gca, 'CLim', cLims);
    if (row == 1)
        xlabel('spatial frequency, X (c/deg)');
    end
    ylabel('spatial frequency, Y (c/deg)');
    title(sprintf('%s (tf = %2.1f Hz)', dataLabel, sData.tfSupport(maxTFindex)));
    
    subplot('Position', subplotPosVectors(1+row,2).v);
    imagesc(sData.tfSupport, sData.spatialFrequencySupport, 10*log10(XTspectra));
    axis 'square';
    set(gca, 'XLim', tfLims, 'YLim', sfLims , 'FontSize', 14);
    set(gca, 'CLim', cLims);
    ylabel('spatial frequency, X (c/deg)');
    if (row == 1)
        xlabel('temporal frequency (Hz)');
    end
    title(sprintf('sfY = %2.1f c/deg', sData.spatialFrequencySupport(maxSFYindex)));
    
    
    subplot('Position', subplotPosVectors(1+row,3).v);
    imagesc(sData.tfSupport, sData.spatialFrequencySupport, 10*log10(YTspectra));
    axis 'square';
    set(gca, 'XLim', tfLims, 'YLim', sfLims , 'FontSize', 14);
    set(gca, 'CLim', cLims);
    ylabel('spatial frequency, Y (c/deg)');
    if (row == 1)
        xlabel('temporal frequency (Hz)');
    end
    title(sprintf('sfX = %2.1f c/deg', sData.spatialFrequencySupport(maxSFXindex)));
    
    subplot('Position', subplotPosVectors(1+row,4).v);
    hold on
    for k = 1:sampledSFsNum 
        targetSF = k-1;
        [~,sfYindex] = min(abs(sData.spatialFrequencySupport-targetSF));
        legends{k} = sprintf('%2.1fcpd', sData.spatialFrequencySupport(sfYindex));
        tfProfile = squeeze(YTspectra(sfYindex,:));
        tfProfile = 0.5*(tfProfile + fliplr(tfProfile));
        plot(sData.tfSupport, squeeze(10*log10(tfProfile)), 'r-', 'LineWidth', 1.5, 'Color', squeeze(lineColors(k,:)));
    end
    axis 'square'
    box on; grid on;
    set(gca, 'XTick', [1 3 10 30 100], 'XScale', 'log', 'XLim', [0.5 max(sData.tfSupport)], 'YLim', cLims, 'FontSize', 14);
    if (row == 1)
        xlabel('temporal frequency (Hz)');
    end
    ylabel('power (dB)');
    hL = legend(legends);
    set(hL, 'units', 'normalized');
    p = get(hL, 'Position');
    p(1) = p(1) + 0.05;
    p(2) = p(2) + 0.01;
    set(hL, 'Position', p);
    
    
    subplot('Position', subplotPosVectors(1+row,5).v);
    hold on
    for k = 1:(sampledTFsNum+1)
        targetHz = 2*(k-1);
        [~,tfIndex] = min(abs(sData.tfSupport-targetHz));
        if (k == sampledTFsNum + 1)
            legends{k} = sprintf('total');
            sfProfile = squeeze(mean(YTspectra,2));
            lineColor = [0 0 0];
        else
            legends{k} = sprintf('%2.1fHz', sData.tfSupport(tfIndex));
            sfProfile = squeeze(YTspectra(:,tfIndex));
            lineColor = squeeze(lineColors(k,:));
        end
        sfProfile = 0.5*(sfProfile + fliplr(sfProfile));
        plot(sData.spatialFrequencySupport, squeeze(10*log10(sfProfile)), 'r-', 'LineWidth', 1.5, 'Color', lineColor);
    end
    axis 'square'
    box on; grid on;
    set(gca, 'XTick', [1 3 10 30 100],'XScale', 'log', 'XLim', [0.5 max(sData.spatialFrequencySupport)], 'YLim', cLims, 'FontSize', 14);
    if (row == 1)
        xlabel('spatial frequency, Y (c/deg)');
    end
    
    hL = legend(legends);
    set(hL, 'units', 'normalized');
    p = get(hL, 'Position');
    p(1) = p(1) + 0.05;
    p(2) = p(2) + 0.01;
    set(hL, 'Position', p);
    
    
    
    subplot('Position', subplotPosVectors(1+row,6).v);
    hold on
    for k = 1:(sampledTFsNum + 1)
        targetHz = 2*(k-1);
        [~,tfIndex] = min(abs(sData.tfSupport-targetHz));
        if (k == sampledTFsNum + 1)
            legends{k} = sprintf('total');
            sfProfile = squeeze(mean(XTspectra,2));
            lineColor = [0 0 0];
        else
            legends{k} = sprintf('%2.1fHz', sData.tfSupport(tfIndex));
            sfProfile = squeeze(XTspectra(:,tfIndex));
            lineColor = squeeze(lineColors(k,:));
        end
        sfProfile = 0.5*(sfProfile + fliplr(sfProfile));
        plot(sData.spatialFrequencySupport, squeeze(10*log10(sfProfile)), 'r-', 'LineWidth', 1.5, 'Color', lineColor);
    end
    axis 'square'
    box on; grid on;
    set(gca, 'XTick', [1 3 10 30 100],'XScale', 'log', 'XLim', [0.5 max(sData.spatialFrequencySupport)], 'YLim', cLims, 'FontSize', 14);
    if (row == 1)
        xlabel('spatial frequency, X (c/deg)');
    end
    hL = legend(legends);
    set(hL, 'units', 'normalized');
    p = get(hL, 'Position');
    p(1) = p(1) + 0.05;
    p(2) = p(2) + 0.01;
    set(hL, 'Position', p);
    
    colormap(hot(1024));
end


