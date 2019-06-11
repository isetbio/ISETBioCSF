function visualizeStabilizedAndDynamicsSpectra(sData, sDataStabilized, figNo)

    % Limits
    sfLims = [0 30];
    tfLims = [-100 100];
    cLims = [-110 45]; % max(sData.meanSpatioTemporalSpectalDensity(:)) + [-40 0];
    
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
    targetYSF = 0;
    [~,sfYIndex] = min(abs(sDataDynamic.spatialFrequencySupport)-targetYSF);
    XTspectraDynamic = squeeze(sDataDynamic.meanSpatioTemporalSpectalDensity(:,sfYIndex,:));
    XTspectraStabilized = squeeze(sDataStabilized.meanSpatioTemporalSpectalDensity(:,sfYIndex,:));
    
    targetSF = 0;
    [~,sfIndex0] = min(abs(sDataDynamic.spatialFrequencySupport)-targetSF);
    
    targetTF = 0;
    [~, tfIndex0] = min(abs(sDataDynamic.tfSupport)-targetTF);
    
    
    sfIndices = (sfIndex0+1):numel(sDataDynamic.spatialFrequencySupport);
    tfIndices = (tfIndex0+1):numel(sDataDynamic.tfSupport);
    
    
    averageXTspectraDynamic = zeros(numel(tfIndices)+1, numel(sfIndices)+1);
    averageXTspectraStabilized = zeros(numel(tfIndices)+1, numel(sfIndices)+1);
    
    averageXTspectraDynamic(1,:) = XTspectraDynamic(tfIndex0,sfIndex0:end);
    averageXTspectraDynamic(:,1) = XTspectraDynamic(tfIndex0:end,sfIndex0);
    averageXTspectraStabilized(1,:) = XTspectraStabilized(tfIndex0,sfIndex0:end);
    averageXTspectraStabilized(:,1) = XTspectraStabilized(tfIndex0:end,sfIndex0);
    
    for sf = 1:numel(sfIndices)
        sfIndex = sfIndices(sf);
        sfIndex2 = sfIndex0 - (sfIndex-sfIndex0);
        for tf = 1:numel(tfIndices)
            tfIndex = tfIndices(tf);
            tfIndex2 = tfIndex0 - (tfIndex-tfIndex0);
            averageXTspectraDynamic(tf+1,sf+1) = 0.25*(...
                XTspectraDynamic(tfIndex, sfIndex) + ...
                XTspectraDynamic(tfIndex2, sfIndex) + ...
                XTspectraDynamic(tfIndex2, sfIndex2) + ...
                XTspectraDynamic(tfIndex, sfIndex2));
            averageXTspectraStabilized(tf+1,sf+1) = 0.25*(...
                XTspectraStabilized(tfIndex, sfIndex) + ...
                XTspectraStabilized(tfIndex2, sfIndex) + ...
                XTspectraStabilized(tfIndex2, sfIndex2) + ...
                XTspectraStabilized(tfIndex, sfIndex2));
        end
    end
    
    
    sfIndices = sfIndex0:numel(sDataDynamic.spatialFrequencySupport);
    tfIndices = tfIndex0:numel(sDataDynamic.tfSupport);
    
    sfAxis = sDataDynamic.spatialFrequencySupport(sfIndices);
    tfAxis = sDataDynamic.tfSupport(tfIndices);
    [X,Y] = meshgrid(sfAxis, tfAxis);
    
    % Do not let 0, so we can plot in log coords
    sfAxis(1) = sfAxis(2)-0.1;
    tfAxis(1) = tfAxis(2)-0.1;

    
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 4, ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.07, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.08, ...
       'topMargin',      0.05);

    sfLims = [sfAxis(2) 80];
    tfLims = [tfAxis(2) 100];
    figure(22); clf;
    subplot('Position', subplotPosVectors(1,1).v);
    [c, h] = contourf(X,Y,10*log10(averageXTspectraDynamic), [cLims(1):2:cLims(end)]);
    set(h,'LineColor','none')
    xlabel('spatial frequency (c/deg)');
    ylabel('temporal frequency (Hz)');
    title('dynamic stimulus');
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', sfLims, 'YLim', tfLims, 'CLim', cLims, 'ZLim', cLims);
    set(gca, 'XTick', [0.1 1 3 10 30 100], 'YTick', [0.1 1 3 10 30 100], 'FontSize', 14);
    axis 'square';
    cBar = colorbar('horizontal');
    cBar.Label.String = 'power (dB)';
    
    subplot('Position', subplotPosVectors(2,1).v);
    tfIndices = find(tfAxis> tfAxis(1));
    [~,tfIndex0] = min(abs(tfAxis));
    plot(sfAxis, 10*log10(squeeze(averageXTspectraDynamic(tfIndex0,:))), 'r-', 'LineWidth', 1.5); hold on
    plot(sfAxis, 10*log10(squeeze(averageXTspectraStabilized(tfIndex0,:))), 'k--', 'LineWidth', 1.5); 
    axis 'square'
    legend({'dynamic', 'stabilized'}, 'Location', 'NorthWest');
    set(gca, 'XScale', 'log',  'XLim', sfLims, 'YLim', [0 45], 'FontSize', 14);
    xlabel('spatial frequency (c/deg)');
    ylabel('power (dB)');
    title('power at tf = 0 Hz');
    
    subplot('Position', subplotPosVectors(1,2).v);
    [c, h] = contourf(X,Y,10*log10(averageXTspectraStabilized), [cLims(1):2:cLims(end)]);
    set(h,'LineColor','none')
    xlabel('spatial frequency (c/deg)');
    ylabel('temporal frequency (Hz)');
    title('stabilized stimulus');
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', sfLims, 'YLim', tfLims, 'CLim', cLims, 'ZLim', cLims);
    set(gca, 'XTick', [0.1 1 3 10 30 100], 'YTick', [0.1 1 3 10 30 100], 'FontSize', 14);
    axis 'square';
    cBar = colorbar('horizontal');
    cBar.Label.String = 'power (dB)';
    

    subplot('Position', subplotPosVectors(2,2).v);
    plot(sfAxis, 10*log10(sum(averageXTspectraDynamic(tfIndices,:),1)), 'r-', 'LineWidth', 1.5); hold on
    plot(sfAxis, 10*log10(sum(averageXTspectraStabilized(tfIndices,:),1)), 'k--', 'LineWidth', 1.5); 
    axis 'square'
    legend({'dynamic', 'stabilized'}, 'Location', 'NorthWest');
    set(gca, 'XScale',  'log', 'XLim', sfLims, 'YLim', [0 45], 'FontSize', 14);
    xlabel('spatial frequency (c/deg)');
    ylabel('power (dB)');
    title('total power at tf > 0 Hz');
    
    sampledTFs = [0 5 10 20 50 100];
    sampledTFsNum = numel(sampledTFs);
    
    sampledSFs = [0 3 10 20 40 60];
    sampledSFsNum = numel(sampledSFs);
    
    lineColors = brewermap(max([sampledSFsNum sampledTFsNum]) , 'Spectral');
    
    legendsTF = {};
    for k = 1:sampledTFsNum
        [~, tfIndex(k)] = min(abs(tfAxis-sampledTFs(k)));
        diffEnergyAsAFunctionOfSF(k,:) =  10*log10(averageXTspectraDynamic(tfIndex(k),:)) - 10*log10(averageXTspectraStabilized(tfIndex(k),:));
        legendsTF{k} = sprintf('%2.0fHz', tfAxis(tfIndex(k)));
    end
    
    legendsSF = {};
    for k = 1:sampledSFsNum
        [~, sfIndex(k)] = min(abs(sfAxis-sampledSFs(k)));
        diffEnergyAsAFunctionOfTF(k,:) =  10*log10(averageXTspectraDynamic(:,sfIndex(k))) - 10*log10(averageXTspectraStabilized(:,sfIndex(k)));
        legendsSF{k} = sprintf('%2.0fc/deg', sfAxis(sfIndex(k)));
    end
    
    subplot('Position', subplotPosVectors(1,3).v); hold on
    for k = 1:sampledTFsNum
        plot(sfAxis, squeeze(diffEnergyAsAFunctionOfSF(k,:)), 'k-', 'Color', squeeze(lineColors(k,:))*0.7, 'LineWidth', 4);
    end
    for k = 1:sampledTFsNum
        plot(sfAxis, squeeze(diffEnergyAsAFunctionOfSF(k,:)), 'k-', 'Color', squeeze(lineColors(k,:)), 'LineWidth', 2);
    end
    axis 'square'
    box on; grid on;
    hL = legend(legendsTF, 'Location', 'NorthWest');
    set(gca, 'XLim', sfLims, 'YLim', [-30 70], 'XScale', 'log', 'FontSize', 14);
    set(gca, 'XTick', [1 3 10 30 60 100]);
    xlabel('spatial frequency (c/deg)');
    ylabel('dynamic-stabilized diff power (dB)');
   
    
    subplot('Position', subplotPosVectors(1,4).v);
    hold on
    for k = 1:sampledSFsNum
        plot(tfAxis, squeeze(diffEnergyAsAFunctionOfTF(k,:)), 'k-', 'Color', squeeze(lineColors(k,:))*0.7, 'LineWidth', 4);
    end
    for k = 1:sampledSFsNum
        plot(tfAxis, squeeze(diffEnergyAsAFunctionOfTF(k,:)), 'k-', 'Color', squeeze(lineColors(k,:)), 'LineWidth', 2);
    end
    axis 'square'
    box on; grid on;
    hL = legend(legendsSF, 'Location', 'NorthWest');
    set(gca, 'XLim', tfLims, 'YLim', [-30 70], 'XScale', 'log', 'FontSize', 14);
    set(gca, 'XTick', [1 3 10 30 100]);
    xlabel('temporal frequency (Hz)');
    ylabel('dynamic-stabilized diff power (dB)');
    colormap(hot)
    
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


