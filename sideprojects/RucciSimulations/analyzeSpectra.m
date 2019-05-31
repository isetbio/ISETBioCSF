function analyzeSpectra(sData, sDataStabilized, figNo)

    % Limits
    sfLims = [0 30];
    tfLims = [-100 100];
    cLims = [-40 30]; % max(sData.meanSpatioTemporalSpectalDensity(:)) + [-40 0];
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1950 550]);
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
    
    
    subplot(2,5,1+row*5);
    imagesc(sData.spatialFrequencySupport, sData.spatialFrequencySupport, 10*log10(XYspectra));
    axis 'square';
    set(gca, 'XLim', max(sfLims)*[-1/2 1/2], 'YLim', sfLims , 'FontSize', 14);
    set(gca, 'CLim', cLims);
    xlabel('spatial frequency, X (c/deg)');
    ylabel('spatial frequency, Y (c/deg)');
    title(sprintf('%s (tf = %2.1f Hz)', dataLabel, sData.tfSupport(maxTFindex)));
    
    subplot(2,5,2+row*5);
    imagesc(sData.tfSupport, sData.spatialFrequencySupport, 10*log10(XTspectra));
    axis 'square';
    set(gca, 'XLim', tfLims, 'YLim', sfLims , 'FontSize', 14);
    set(gca, 'CLim', cLims);
    ylabel('spatial frequency, X (c/deg)');
    xlabel('temporal frequency (Hz)');
    title(sprintf('sfY = %2.1f c/deg', sData.spatialFrequencySupport(maxSFYindex)));
    
    
    subplot(2,5,3+row*5)
    imagesc(sData.tfSupport, sData.spatialFrequencySupport, 10*log10(YTspectra));
    axis 'square';
    set(gca, 'XLim', tfLims, 'YLim', sfLims , 'FontSize', 14);
    set(gca, 'CLim', cLims);
    ylabel('spatial frequency, Y (c/deg)');
    xlabel('temporal frequency (Hz)');
    title(sprintf('sfX = %2.1f c/deg', sData.spatialFrequencySupport(maxSFXindex)));
    
    subplot(2,5,4+row*5);
    hold on
    for k = 1:sampledSFsNum 
        targetSF = k-1;
        [~,sfYindex] = min(abs(sData.spatialFrequencySupport-targetSF));
        legends{k} = sprintf('%2.1f c/deg', sData.spatialFrequencySupport(sfYindex));
        tfProfile = squeeze(YTspectra(sfYindex,:));
        tfProfile = 0.5*(tfProfile + fliplr(tfProfile));
        plot(sData.tfSupport, squeeze(10*log10(tfProfile)), 'r-', 'LineWidth', 1.5, 'Color', squeeze(lineColors(k,:)));
    end
    axis 'square'
    box on; grid on;
        hL = legend(legends, 'Location', 'eastoutside');
        %set(hL,'units','normalized');
        %set(hL,'position',[0.47,0.1,0.1,0.18]);
     
    set(gca, 'XTick', [1 3 10 30 100], 'XScale', 'log', 'XLim', [0.5 max(sData.tfSupport)], 'YLim', cLims, 'FontSize', 14);
    xlabel('temporal frequency (Hz)');
    ylabel('power (dB)');
    
    subplot(2,5,5+row*5);
    hold on
    for k = 1:sampledTFsNum 
        targetHz = 2*(k-1);
        [~,tfIndex] = min(abs(sData.tfSupport-targetHz));
        legends{k} = sprintf('%2.1f Hz', sData.tfSupport(tfIndex));
        sfProfile = squeeze(YTspectra(:,tfIndex));
        sfProfile = 0.5*(sfProfile + fliplr(sfProfile));
        plot(sData.spatialFrequencySupport, squeeze(10*log10(sfProfile)), 'r-', 'LineWidth', 1.5, 'Color', squeeze(lineColors(k,:)));
    end
    axis 'square'
    box on; grid on;
        hL = legend(legends, 'Location', 'eastoutside');
        %set(hL,'units','normalized');
        %set(hL,'position',[0.47,0.1,0.1,0.18]);
    
    set(gca, 'XTick', [1 3 10 30 100],'XScale', 'log', 'XLim', [0.5 max(sData.spatialFrequencySupport)], 'YLim', cLims, 'FontSize', 14);
    xlabel('spatial frequency, Y (c/deg)');
    ylabel('power (dB)');
    
    colormap(hot(1024));
end


