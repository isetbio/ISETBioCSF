function visualizeStabilizedAndDynamicsSpectra(sData, sDataStabilized, figNo)

    % Limits
    sfLims = [0 30];
    tfLims = [-100 100];
    dbRange = [-20 65]; %  [-60 45]; % in dB
    
    hFig = figure(figNo+100); clf;
    plotSummarySlices(sDataStabilized, sData, dbRange );
    pause
    
    hFig = figure(figNo+100); clf;
    set(hFig, 'Position', [10 10 2250 625]);
    plotSlices(sDataStabilized,'STABILIZED', sfLims, tfLims, dbRange ,0);
    plotSlices(sData,'DYNAMIC', sfLims, tfLims, dbRange ,1);
    drawnow;
end

function plotSummarySlices(sDataStabilized, sDataDynamic,  dbRange )

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
    
    % Report total power
    totalPowerOfDynamicStimulusDB = 10*log10(sum(averageXTspectraDynamic(:)))
    totalPowerOfStabilizedStimulusDB = 10*log10(sum(averageXTspectraStabilized(:)))
    
    % Compute sums along SF, TF
    sumStabilizedAcrossSFs = sum(averageXTspectraStabilized, 2);
    sumStabilizedAcrossTFs = sum(averageXTspectraStabilized, 1);
    sumDynamicAcrossSFs = sum(averageXTspectraDynamic, 2);
    sumDynamicAcrossTFs = sum(averageXTspectraDynamic, 1);
     
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 4, ...
       'heightMargin',  0.1, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.04, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.08, ...
       'topMargin',      0.03);

    logPlot = true;
    if (logPlot)
        % make the zero frequency point 0.1 so we can display it on log plot
        sfSupport(1) = 0.1;
        tfSupport(1) = 0.1;
    end
    
    sfLims = [sfSupport(1) 50];
    tfLims = [tfSupport(1) 100];
    
    
    figure(22); clf;
    
    % The spatiotemporal spectrum of the stabilized stimulus
    ax = subplot('Position', subplotPosVectors(1,1).v);
    visualize2Dspectrum(ax, sfSupport, tfSupport, averageXTspectraStabilized, sfLims, tfLims, dbRange, 'stabilized', logPlot);

    
    % The spatiotemporal spectrum of the dynamic stimulus
    ax = subplot('Position', subplotPosVectors(1,2).v);
    visualize2Dspectrum(ax, sfSupport, tfSupport, averageXTspectraDynamic, sfLims, tfLims, dbRange, 'dynamic', logPlot);

    
    
    % Visualize summed spectrum across the SF axis 
    subplot('Position', subplotPosVectors(1,3).v); 
    visualizeMarginalSpectra(tfSupport, sumStabilizedAcrossSFs, sumDynamicAcrossSFs, tfLims, 'temporal frequency (Hz)', logPlot);
    
    
    % Visualized summed spectrum across the TF axis
    subplot('Position', subplotPosVectors(1,4).v); 
    visualizeMarginalSpectra(sfSupport, sumStabilizedAcrossTFs, sumDynamicAcrossTFs, sfLims, 'spatial frequency (c/deg)', logPlot);
    
    
    sampledTFs = [0 2 5 7 10 20 50 100];
    sampledTFsNum = numel(sampledTFs);
    
    sampledSFs = [0 3 10 20 40 60];
    sampledSFsNum = numel(sampledSFs);
    lineColors = brewermap(max([sampledSFsNum sampledTFsNum]) , 'Spectral');
    
    
    % Spectral slices at different TFs of the stabilized stimulus
    subplot('Position', subplotPosVectors(2,1).v);
    samplingDimension = 1;
    visualizeSpectralSlices(sfSupport, sfLims, 'spatial frequency, X (c/deg)', ...
        tfSupport, sampledTFs, 'Hz', samplingDimension, averageXTspectraStabilized, dbRange, lineColors, 'stabilized', logPlot);
    
    
    % Spectral (sf) slices at different TFs of the dynamic stimulus
    subplot('Position', subplotPosVectors(2,2).v);
    visualizeSpectralSlices(sfSupport, sfLims, 'spatial frequency, X (c/deg)', ...
        tfSupport, sampledTFs, 'Hz', samplingDimension, averageXTspectraDynamic, dbRange, lineColors, 'dynamic', logPlot);

    
%     % Spectral (sf) slices of differential power (dynamic-stabilized) stimulus
%     subplot('Position', subplotPosVectors(2,3).v); hold on
%     % integrate over TF
%     tfBands = find(tfSupport >= sfSupport(1) & tfSupport < 100);
%     totalDynamicPower = sum(averageXTspectraDynamic(tfBands,:),1);
%     totalStabilizedPower = sum(averageXTspectraStabilized,1);
%     plot(sfSupport, 10*log10(totalDynamicPower), 'r-', 'LineWidth', 1.5);
%     hold on
%     plot(sfSupport, 10*log10(totalStabilizedPower), 'k-', 'LineWidth', 1.5);
%     hL = legend({'dynamic', 'stabilized'},  'Location', 'northeast');
%     set(gca, 'XLim', sfLims, 'YLim', cLims, 'FontSize', 14);
%     if (logPlot)
%         set(gca, 'XScale', 'log');
%     end
%     set(gca, 'XTick', [0.1 1 3 10 30 60 100]);
%     xlabel('spatial frequency, X (c/deg)');
%     ylabel('power (dB)');
%     box on; grid on
    
end

function [averageSpectra, sfAxis, tfAxis] = averageAcrossQuadrantsSpectum(sfSupport, tfSupport, xtSpectra)

    indexOfZeroSF = floor(size(xtSpectra,2)/2)+1;
    indexOfZeroTF = floor(size(xtSpectra,1)/2)+1;
    
    sfIndices = indexOfZeroSF:size(xtSpectra,2);
    tfIndices = indexOfZeroTF:size(xtSpectra,1);
    averageSpectra = zeros(numel(tfIndices), numel(sfIndices));
    
    sfAxis = sfSupport(sfIndices);
    tfAxis = tfSupport(tfIndices);

    for tf = 1:numel(tfIndices)
        tfIndex = tfIndices(tf);
        tfIndexSymmetric = indexOfZeroTF - (tfIndex - indexOfZeroTF);
        for sf = 1:numel(sfIndices)
            sfIndex = sfIndices(sf);
            sfIndexSymmetric = indexOfZeroSF - (sfIndex - indexOfZeroSF);
            if ((tf == 1) && (sf == 1))
                % (0,0) frequency
                averageSpectra(tf,sf) = xtSpectra(tfIndex ,sfIndex);
            else
                if (tf == 1) && (sf > 1)
                    % at zero TF axis, sum corresponding points from 2
                    % spatial frequencies
                    averageSpectra(tf,sf) = ...
                        xtSpectra(tfIndex,sfIndex) + ...
                        xtSpectra(tfIndex,sfIndexSymmetric);
                elseif (sf == 1) && (tf > 1)
                    % at zero SF axis, sum corresponding points from 2
                    % temporal frequencies
                    averageSpectra(tf,sf) = ...
                        xtSpectra(tfIndex,sfIndex) + ...
                        xtSpectra(tfIndexSymmetric,sfIndex);
                else
                    % sum corresponding points from 4 quadrants
                    averageSpectra(tf,sf) = ...
                        xtSpectra(tfIndex,sfIndex) + ...
                        xtSpectra(tfIndex,sfIndexSymmetric) + ...
                        xtSpectra(tfIndexSymmetric,sfIndex) + ...
                        xtSpectra(tfIndexSymmetric,sfIndexSymmetric);
                end
            end
        end
    end 
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


function visualize2Dspectrum(ax, sfSupport, tfSupport, spectrum2D, sfLims, tfLims, dbRange, figureTitle, logPlot )
    [X,Y] = meshgrid(sfSupport, tfSupport);
    spectrum2D = 10*log10(spectrum2D);
    spectrum2D(spectrum2D<dbRange(1)) = dbRange(1);
    spectrum2D(spectrum2D>dbRange(2)) = dbRange(2);
    dbLevels = linspace(dbRange(1), dbRange(2),100);
    [~, h] = contourf(X,Y,spectrum2D, dbLevels);
    set(h,'LineColor','none')
    xlabel('\it spatial frequency, X (c/deg)');
    ylabel('\it temporal frequency (Hz)');
    title(figureTitle);
    set(gca, 'XLim', sfLims, 'YLim', tfLims, 'CLim', dbRange, 'ZLim', dbRange);
    if (logPlot)
        set(gca, 'XScale', 'log', 'YScale', 'log');
    end
    set(gca, 'XTick', [0.1 1 3 10 30 100], 'YTick', [0.1 1 3 10 30 100], 'FontSize', 14);
    cBar = colorbar('horizontal', 'Location', 'northoutside');
    cBar.Label.String = 'power (dB)';
    colormap(ax, hot(1024));
    drawnow;
end

function visualizeSpectralSlices(xAxisSupport, xAxisLims, xAxisLabel, sampledAxisSupport, sampledPoints, samplePointUnit, samplingDimension, spectrum2D, dbRange, ...
    lineColors, figureTitle, logPlot)
    legends = {};
    for k = 1:numel(sampledPoints)
        [~,idx] = min(abs(sampledAxisSupport-sampledPoints(k)));
        if (samplingDimension == 1)
            slice = 10*log10(spectrum2D(idx,:)); 
        else
            slice = 10*log10(spectrum2D(:,idx));
        end
        if (numel(xAxisSupport) ~= numel(slice))
            error('inconsistent sampling dimension');
        end
        plot(xAxisSupport, slice, 'k-', 'Color', squeeze(lineColors(k,:)), 'LineWidth', 2);
        hold on;
        legends{k} = sprintf('%2.2f %s', sampledAxisSupport(idx), samplePointUnit);
    end
    
    %axis 'square'
    box on; grid on;
    hL = legend(legends, 'NumColumns',2, 'Location', 'northeast');
    set(gca, 'XLim', xAxisLims, 'YLim', dbRange, 'FontSize', 14);
    if (logPlot)
        set(gca, 'XScale', 'log');
    end
    set(gca, 'XTick', [1 3 10 30 60 100]);
    xlabel(xAxisLabel);
    ylabel('power (dB)');
    title(figureTitle);
end

function visualizeMarginalSpectra(support, sumStabilized, sumDynamic, xAxisLims, xAxisLabel, logPlot)
    plot(support, 10*log10(sumStabilized), 'k-', 'LineWidth', 1.5);
    hold on;
    plot(support, 10*log10(sumDynamic), 'r--', 'LineWidth', 1.5);
    plot(support, 10*log10(sumDynamic)-10*log10(sumStabilized), 'b-', 'LineWidth', 1.5);
    hL = legend({'stabilized', 'dynamic', 'dynamic-stabilized'},  'Location', 'northeast');
    set(gca, 'XLim', xAxisLims,'FontSize', 14);
    if (logPlot)
        set(gca, 'XScale', 'log');
    end
    set(gca, 'XTick', [0 1 3 10 30 60 100]);
    xlabel(xAxisLabel);
    ylabel('power (dB)');
    box on; grid on
end
    