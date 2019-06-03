function analyzeSpectra(sData, sDataStabilized, figNo)

    % Limits
    sfLims = [0 30];
    tfLims = [-100 100];
    cLims = [-60 35]; % max(sData.meanSpatioTemporalSpectalDensity(:)) + [-40 0];
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 2250 625]);
    plotSlices(sDataStabilized,'STABILIZED', sfLims, tfLims, cLims,0);
    plotSlices(sData,'DYNAMIC', sfLims, tfLims, cLims,1);
    
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


