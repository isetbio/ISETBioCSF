function visualizeStabilizedAndDynamicsSpectra(sData, sDataStabilized, figNo)

    % Limits
    sfLims = [0 50];
    tfLims = [0 200];
    dbRange = [-20 65]; %  [-60 45]; % in dB
    
    hFig = figure(figNo+100); clf;
    plotSummarySlices(sDataStabilized, sData, dbRange, sfLims, tfLims);
end

function plotSummarySlices(sDataStabilized, sDataDynamic,  dbRange, sfLims, tfLims)

    % organization of sDataStabilized.meanSpatioTemporalPowerSpectalDensity is [TF, SFy, SFx]
    
    spatiotemporalSpectrumMethod = 'select sfY=0';
    spatiotemporalSpectrumMethod = 'sum over sfY';
    %spatiotemporalSpectrumMethod = 'radial averaging';
    
    switch (spatiotemporalSpectrumMethod)
        case 'sum over sfY'
            % Sum over the ySF axis
            XTpowerSpectralDensityStabilized = squeeze(sum(sDataStabilized.meanSpatioTemporalPowerSpectalDensity,2));
            XTpowerSpectralDensityDynamic = squeeze(sum(sDataDynamic.meanSpatioTemporalPowerSpectalDensity,2));
            
            % Average along 4 quadrants
            [averageXTpowerSpectralDensityDynamic, sfSupport, tfSupport] = ...
                averageAcrossQuadrantsSpectum(sDataDynamic.spatialFrequencySupport, sDataDynamic.tfSupport, XTpowerSpectralDensityDynamic);
            [averageXTpowerSpectralDensityStabilized, ~,~] = ...
                averageAcrossQuadrantsSpectum(sDataDynamic.spatialFrequencySupport, sDataDynamic.tfSupport, XTpowerSpectralDensityStabilized);
            
        case 'select sfY=0'
            % Select the ySF = 0 slice
            targetSF = 0;
            [~,sfYIndex] = min(abs(sDataDynamic.spatialFrequencySupport-targetSF));
            XTpowerSpectralDensityStabilized = squeeze(sDataStabilized.meanSpatioTemporalPowerSpectalDensity(:,sfYIndex,:));
            XTpowerSpectralDensityDynamic = squeeze(sDataDynamic.meanSpatioTemporalPowerSpectalDensity(:,sfYIndex,:));
            
            % Average along 4 quadrants
            [averageXTpowerSpectralDensityDynamic, sfSupport, tfSupport] = ...
                averageAcrossQuadrantsSpectum(sDataDynamic.spatialFrequencySupport, sDataDynamic.tfSupport, XTpowerSpectralDensityDynamic);
            [averageXTpowerSpectralDensityStabilized, ~,~] = ...
                averageAcrossQuadrantsSpectum(sDataDynamic.spatialFrequencySupport, sDataDynamic.tfSupport, XTpowerSpectralDensityStabilized);
    
        case 'radial averaging' 
            m = numel(sDataDynamic.spatialFrequencySupport)/2+1;
            sfSupport = sDataDynamic.spatialFrequencySupport(m:end);
            m = numel(sDataDynamic.tfSupport)/2+1;
            tfSupport = sDataDynamic.tfSupport(m:end);
    
            indexOfZeroTF = floor(size(sDataDynamic.meanSpatioTemporalPowerSpectalDensity,1)/2)+1;
            tfIndices = indexOfZeroTF:size(sDataDynamic.meanSpatioTemporalPowerSpectalDensity,1);

            for tf = 1:numel(tfIndices)
                tfIndex = tfIndices(tf);
                tfIndexSymmetric = indexOfZeroTF - (tfIndex - indexOfZeroTF);
                if (tf == 1)
                    spectrumXY = squeeze(sDataDynamic.meanSpatioTemporalPowerSpectalDensity(tfIndex,:,:));
                else
                    spectrumXY = squeeze(sDataDynamic.meanSpatioTemporalPowerSpectalDensity(tfIndex,:,:)) + ...
                             squeeze(sDataDynamic.meanSpatioTemporalPowerSpectalDensity(tfIndexSymmetric,:,:));
                end
                meanSlice = radiallyAveragedSpectrum(spectrumXY);
                
                if (tf == 1)
                    averageXTpowerSpectralDensityDynamic = zeros(numel(tfIndices), length(meanSlice));
                    averageXTpowerSpectralDensityStabilized = zeros(numel(tfIndices), length(meanSlice));
                end
                averageXTpowerSpectralDensityDynamic(tf,:) = meanSlice;
                if (tf == 1)
                    spectrumXY = squeeze(sDataStabilized.meanSpatioTemporalPowerSpectalDensity(tfIndex,:,:));
                else
                    spectrumXY = squeeze(sDataStabilized.meanSpatioTemporalPowerSpectalDensity(tfIndex,:,:)) + ...
                                squeeze(sDataStabilized.meanSpatioTemporalPowerSpectalDensity(tfIndexSymmetric,:,:));
                end
                averageXTpowerSpectralDensityStabilized(tf,:) = radiallyAveragedSpectrum(spectrumXY);
            end
        otherwise
            error('Unknown XT spectral method: ''%s''.', spatiotemporalSpectrumMethod)
    end

    
    % Report total power
    deltaSF = sfSupport(2)-sfSupport(1);
    deltaTF = tfSupport(2)-tfSupport(1);
    totalPowerOfDynamicStimulusDB = 10*log10(deltaSF*deltaTF*sum(averageXTpowerSpectralDensityDynamic(:)))
    totalPowerOfStabilizedStimulusDB = 10*log10(deltaSF*deltaTF*sum(averageXTpowerSpectralDensityStabilized(:)))
    pause
    
    % Compute sums along SF, TF
    sumStabilizedAcrossSFs = sum(averageXTpowerSpectralDensityStabilized, 2);
    sumStabilizedAcrossTFs = sum(averageXTpowerSpectralDensityStabilized, 1);
    sumDynamicAcrossSFs = sum(averageXTpowerSpectralDensityDynamic, 2);
    sumDynamicAcrossTFs = sum(averageXTpowerSpectralDensityDynamic, 1);
     
    idx = find(tfSupport < 0.4);
    sumDynamicAcrossSFsLowTF = sum(sumDynamicAcrossSFs(idx));
    sumStabilizedAcrossSFsLowTF = sum(sumStabilizedAcrossSFs(idx));
    diffLowTF = sumDynamicAcrossSFsLowTF-sumStabilizedAcrossSFsLowTF
    
    idx = find(tfSupport >=0.4);
    sumDynamicAcrossSFsHighTF = sum(sumDynamicAcrossSFs(idx));
    sumStabilizedAcrossSFsHighTF = sum(sumStabilizedAcrossSFs(idx));
    diffHighTF = sumDynamicAcrossSFsHighTF - sumStabilizedAcrossSFsHighTF
    
    [sum(sumDynamicAcrossSFs) sum(sumStabilizedAcrossSFs)]
    
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
    
    sfLims = [sfSupport(1) sfLims(2)];
    tfLims = [tfSupport(1) tfLims(2)];
    
    
    figure(22); clf;
    
    % The spatiotemporal spectrum of the stabilized stimulus
    ax = subplot('Position', subplotPosVectors(1,1).v);
    visualize2Dspectrum(ax, sfSupport, tfSupport, averageXTpowerSpectralDensityStabilized, sfLims, tfLims, dbRange, 'stabilized', logPlot);

    
    % The spatiotemporal spectrum of the dynamic stimulus
    ax = subplot('Position', subplotPosVectors(1,2).v);
    visualize2Dspectrum(ax, sfSupport, tfSupport, averageXTpowerSpectralDensityDynamic, sfLims, tfLims, dbRange, 'dynamic', logPlot);

    
     % Visualized summed spectrum across the TF axis
    subplot('Position', subplotPosVectors(1,3).v); 
    marginalLimsDB = [-20 75];
    visualizeMarginalSpectra(sfSupport, sumStabilizedAcrossTFs, sumDynamicAcrossTFs, sfLims, marginalLimsDB, [tfSupport(1) tfSupport(end)], 'Hz', 'spatial frequency (c/deg)', logPlot);
    
    
    % Visualize summed spectrum across the SF axis 
    subplot('Position', subplotPosVectors(1,4).v); 
    visualizeMarginalSpectra(tfSupport, sumStabilizedAcrossSFs, sumDynamicAcrossSFs, tfLims, marginalLimsDB, [sfSupport(1) sfSupport(end)], 'c/deg', 'temporal frequency (Hz)', logPlot);
    
    
    
    sampledTFs = [0 2 5 7 10 20 50 100];
    sampledTFsNum = numel(sampledTFs);
    
    sampledSFs = [0 3 10 20 40 60];
    sampledSFsNum = numel(sampledSFs);
    lineColors = brewermap(max([sampledSFsNum sampledTFsNum]) , 'Spectral');
    
    
    % Spectral slices at different TFs of the stabilized stimulus
    subplot('Position', subplotPosVectors(2,1).v);
    samplingDimension = 1;
    visualizeSpectralSlices(sfSupport, sfLims, 'spatial frequency, X (c/deg)', ...
        tfSupport, sampledTFs, 'Hz', samplingDimension, averageXTpowerSpectralDensityStabilized, dbRange, lineColors, 'stabilized', logPlot);
    
    
    % Spectral (sf) slices at different TFs of the dynamic stimulus
    subplot('Position', subplotPosVectors(2,2).v);
    visualizeSpectralSlices(sfSupport, sfLims, 'spatial frequency, X (c/deg)', ...
        tfSupport, sampledTFs, 'Hz', samplingDimension, averageXTpowerSpectralDensityDynamic, dbRange, lineColors, 'dynamic', logPlot);

    subplot('Position', subplotPosVectors(2,3).v);
    tfRangeOverWhichToSumSlices = [1 30];
    legends = {};
    legends = visualizeSummedSpectralSlices(sfSupport, sfLims, 'spatial frequency, X (c/deg)', ...
        tfSupport, tfRangeOverWhichToSumSlices, 'Hz', samplingDimension, averageXTpowerSpectralDensityStabilized, dbRange,  [0 0 0], legends, 'stabilized', logPlot);
    
    %subplot('Position', subplotPosVectors(2,4).v);
    legends = visualizeSummedSpectralSlices(sfSupport, sfLims, 'spatial frequency, X (c/deg)', ...
        tfSupport, tfRangeOverWhichToSumSlices, 'Hz', samplingDimension, averageXTpowerSpectralDensityDynamic, dbRange, [1 0 0], legends,'dynamic', logPlot);
    
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

    for k = 1:numel(sampledPoints)
        [~,idx] = min(abs(sampledAxisSupport-sampledPoints(k)));
        if (samplingDimension == 1)
            slices(k,:) = 10*log10(spectrum2D(idx,:)); 
        else
            slices(k,:) = 10*log10(spectrum2D(:,idx));
        end
    end
    
    visualizeSlices(xAxisSupport, xAxisLims, xAxisLabel, sampledAxisSupport, sampledPoints, samplePointUnit, slices, dbRange, ...
        lineColors, figureTitle, logPlot);
end

function legends = visualizeSummedSpectralSlices(xAxisSupport, xAxisLims, xAxisLabel, sampledAxisSupport, sampledPointsRange, samplePointUnit, samplingDimension, spectrum2D, dbRange, ...
    lineColor, legends, setTitle, logPlot)

    theZeroSlice = squeeze(spectrum2D(1,:));
    indicesOfSampledSlices = find((sampledAxisSupport >= sampledPointsRange(1)) & (sampledAxisSupport <= sampledPointsRange(2)));
    theSummedSlices = sum(spectrum2D(indicesOfSampledSlices,:),samplingDimension);
    
    plot(xAxisSupport, 10*log10(theZeroSlice), 'k--', 'LineWidth', 3, 'Color', lineColor*0.8); hold on
    plot(xAxisSupport, 10*log10(theSummedSlices), 'k-', 'LineWidth', 3, 'Color', lineColor*0.8);
    legends{numel(legends)+1} = sprintf('%2.1f%s (%s)', sampledAxisSupport(1), samplePointUnit, setTitle);
    legends{numel(legends)+1} = sprintf('%2.1f - %2.1f %s (%s)', sampledAxisSupport(indicesOfSampledSlices(1)), sampledAxisSupport(indicesOfSampledSlices(end)), samplePointUnit, setTitle);
    
    %axis 'square'
    box on; grid on;
    hL = legend(legends, 'Location', 'northeast');
    set(gca, 'XLim', xAxisLims, 'YLim', dbRange, 'FontSize', 14);
    if (logPlot)
        set(gca, 'XScale', 'log');
    end
    set(gca, 'XTick', [1 3 10 30 60 100]);
    xlabel(xAxisLabel);
    ylabel('power (dB)');
    
end


function visualizeDifferentialSpectralSlices(xAxisSupport, xAxisLims, xAxisLabel, sampledAxisSupport, sampledPoints, samplePointUnit, samplingDimension, spectrum2DA, spectrum2DB, dbRange, ...
    lineColors, figureTitle, logPlot)

    for k = 1:numel(sampledPoints)
        [~,idx] = min(abs(sampledAxisSupport-sampledPoints(k)));
        if (samplingDimension == 1)
            slices(k,:) = 10*log10(spectrum2DA(idx,:)) - 10*log10(spectrum2DB(idx,:)); 
        else
            slices(k,:) = 10*log10(spectrum2DA(:,idx)) - 10*log10(spectrum2DB(:,idx));
        end
    end
    
    visualizeSlices(xAxisSupport, xAxisLims, xAxisLabel, sampledAxisSupport, sampledPoints, samplePointUnit, slices, dbRange, ...
        lineColors, figureTitle, logPlot);
    
end

function visualizeSlices(xAxisSupport, xAxisLims, xAxisLabel, sampledAxisSupport, sampledPoints, samplePointUnit,  slices, dbRange, ...
    lineColors, figureTitle, logPlot)
    legends = {};
    for k = 1:numel(sampledPoints)
        [~,idx] = min(abs(sampledAxisSupport-sampledPoints(k)));
        slice = squeeze(slices(k,:));
        if (numel(xAxisSupport) ~= numel(slice))
            error('inconsistent sampling dimension');
        end
        plot(xAxisSupport, slice, 'k-', 'Color', squeeze(lineColors(k,:))*0.8, 'LineWidth', 3);
        hold on;
        legends{k} = sprintf('%2.1f %s', sampledAxisSupport(idx), samplePointUnit);
    end
    
    for k = 1:numel(sampledPoints)
        [~,idx] = min(abs(sampledAxisSupport-sampledPoints(k)));
        slice = squeeze(slices(k,:));
        if (numel(xAxisSupport) ~= numel(slice))
            error('inconsistent sampling dimension');
        end
        plot(xAxisSupport, slice, 'k-', 'Color', squeeze(lineColors(k,:)), 'LineWidth', 1.5);
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

function visualizeMarginalSpectra(support, sumStabilized, sumDynamic, xAxisLims, dbRange, integrationLimits, integrationAxisUnit, xAxisLabel, logPlot)
    hold on;
    plot(support, 10*log10(sumStabilized), 'k-', 'LineWidth',3, 'Color', [ 0 0 0]);
    plot(support, 10*log10(sumDynamic), 'r-', 'LineWidth', 3);
    plot(support, 10*log10(sumDynamic)-10*log10(sumStabilized), 'c-', 'LineWidth', 3);
    plot(support, 10*log10(sumDynamic)-10*log10(sumStabilized), 'b-', 'LineWidth', 3);
    plot(support, 10*log10(sumDynamic)-10*log10(sumStabilized), 'c-', 'LineWidth', 1.5);
    hL = legend({'stabilized', 'dynamic', 'dynamic-stabilized'},  'Location', 'southwest');
    set(gca, 'XLim', xAxisLims, 'YLim', dbRange, 'FontSize', 14);
    if (logPlot)
        set(gca, 'XScale', 'log');
    end
    set(gca, 'XTick', [0 1 3 10 30 60 100]);
    xlabel(xAxisLabel);
    ylabel('power (dB)');
    title(sprintf('integrated power from %2.1f-%2.1f %s', integrationLimits(1), integrationLimits(2), integrationAxisUnit));
    box on; grid on
end
    