function analyzeSpectra(sData, sDataStabilized, figNo)

    % Limits
    sfLims = [0 15];
    tfLims = [-50 50];
    cLims = [-30 20]; % max(sData.meanSpatioTemporalSpectalDensity(:)) + [-40 0];
    
    figure(figNo); clf;
    plotSlices(sDataStabilized,'STABILIZED', sfLims, tfLims, cLims,0);
    plotSlices(sData,'DYNAMIC', sfLims, tfLims, cLims,1);
    
end

function plotSlices(sData, dataLabel, sfLims, tfLims, cLims, column)

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
    lineColors = brewermap(sampledSFsNum , 'Spectral');
    
    
    subplot(4,2,1+column);
    imagesc(sData.spatialFrequencySupport, sData.spatialFrequencySupport, 10*log10(XYspectra));
    axis 'square';
    set(gca, 'XLim', max(sfLims)*[-1/2 1/2], 'YLim', sfLims , 'FontSize', 14);
    set(gca, 'CLim', cLims);
    xlabel('spatial frequency, X (c/deg)');
    ylabel('spatial frequency, Y (c/deg)');
    title(sprintf('%s (tf = %2.1f Hz)', dataLabel, sData.tfSupport(maxTFindex)));
    
    subplot(4,2,3+column);
    imagesc(sData.tfSupport, sData.spatialFrequencySupport, 10*log10(XTspectra));
    axis 'square';
    set(gca, 'XLim', tfLims, 'YLim', sfLims , 'FontSize', 14);
    set(gca, 'CLim', cLims);
    ylabel('spatial frequency, X (c/deg)');
    xlabel('temporal frequency (Hz)');
    title(sprintf('sfY = %2.1f c/deg', sData.spatialFrequencySupport(maxSFYindex)));
    
    
    subplot(4,2,5+column)
    imagesc(sData.tfSupport, sData.spatialFrequencySupport, 10*log10(YTspectra));
    axis 'square';
    set(gca, 'XLim', tfLims, 'YLim', sfLims , 'FontSize', 14);
    set(gca, 'CLim', cLims);
    ylabel('spatial frequency, Y (c/deg)');
    xlabel('temporal frequency (Hz)');
    title(sprintf('sfX = %2.1f c/deg', sData.spatialFrequencySupport(maxSFXindex)));
    
    subplot(4,2,7+column);
    hold on
    for k = 1:sampledSFsNum 
        [~,sfYindex] = min(abs(sData.spatialFrequencySupport-(k-1)));
        legends{k} = sprintf('%2.1f c/deg', sData.spatialFrequencySupport(sfYindex));
        plot(sData.tfSupport, squeeze(10*log10(YTspectra(sfYindex,:))), 'r-', 'LineWidth', 1.5, 'Color', squeeze(lineColors(k,:)));
    end
    box on; grid on;
    hL = legend(legends);
    set(hL,'units','normalized');
    set(hL,'position',[0.47,0.1,0.1,0.18]);
    set(gca, 'XLim', tfLims, 'YLim', cLims, 'FontSize', 14);
    xlabel('temporal frequency (Hz)');
    ylabel('power (dB)');
    
%     subplot(4,2,7+column);
%     powerIntegratedOverTF = squeeze(sum(YTspectra,2));
%     plot(sData.spatialFrequencySupport, 10*log10(powerIntegratedOverTF), 'k-');
%     set(gca, 'XLim', sfLims, 'FontSize', 14);
%     xlabel('spatial frequency, Y (c/deg)');
%     ylabel('power (dB)');
    colormap(hot(1024));
end


