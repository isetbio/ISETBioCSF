function testPhotocurrentNoise

    % Generate desired spectral power distribution
    sampleTime = 0.5/1000;
    totalTimeSeconds = 0.15*4;
    nInstances = 1024;
    
    nSamples = round(totalTimeSeconds/sampleTime);
    [freqAxis, noiseSPD, noiseSPDlowFreq,  ...
        noiseSPDhighFreq, cornerFreqLow, cornerFreqHigh] = ...
        generateSPD(nSamples, sampleTime);

    % Generate noise with desired spectrum
    [timeAxis, noise, realizedNoiseSpectrum] = ...
        generateNoise(nSamples, sampleTime, noiseSPD, nInstances);

    freqRange = [1 500];
    noiseSPDRange = [0.001 0.25];
    visualizeTotalTimeSeconds = 150/1000;
    visualizeAmplitudeRange = 12;
    
    displayResults(timeAxis, noise, freqAxis, noiseSPD, noiseSPDlowFreq, noiseSPDhighFreq, ...
        cornerFreqLow, cornerFreqHigh, realizedNoiseSpectrum, ...
        freqRange, noiseSPDRange, visualizeTotalTimeSeconds, visualizeAmplitudeRange);
end

function  displayResults(timeAxis, noise, freqAxis, noiseSPD, noiseSPDlowFreq, noiseSPDhighFreq, ...
        cornerFreqLow, cornerFreqHigh, realizedNoiseSpectrum, ...
        freqRange, noiseSPDRange, visualizeTotalTimeSeconds, visualizeAmplitudeRange)
    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 650 320], 'Color', [1 1 1]);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', 2, ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.07, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.15, ...
       'topMargin',      0.1);
    
    subplot('Position', subplotPosVectors(1,1).v);
    plot(freqAxis, noiseSPDlowFreq, 'r-', 'LineWidth', 1.5); hold on;
    plot(freqAxis, noiseSPDhighFreq, 'b-', 'LineWidth', 1.5);
    plot(freqAxis, noiseSPD, 'm-', 'LineWidth', 1.5);
    plot(freqAxis, realizedNoiseSpectrum(1:numel(freqAxis)), 'k-', 'LineWidth', 1.5);
    plot(cornerFreqLow, noiseSPDlowFreq(1)/2, 'ro', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5]);
    plot(cornerFreqHigh, noiseSPDhighFreq(1)/2, 'bo', 'LineWidth', 1.5, 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 1]);
    set(gca, 'XScale' ,'log', 'YScale', 'log', 'XLim', freqRange, 'YLim', noiseSPDRange, ...
        'XTick', [0.3 1 3 10 30 100 300]);
    xlabel('\it frequency (Hz)');
    ylabel('\it power spectrum (pAmps^2)');
    set(gca, 'FontSize', 14);
    axis 'square'; grid on; box on;
    
    subplot('Position', subplotPosVectors(1,2).v);
    sampleTime = timeAxis(2)-timeAxis(1);
    N = round(visualizeTotalTimeSeconds/sampleTime);
    timeAxis = timeAxis(1:N)*1000;
    hPlot = plot(timeAxis(1:N) - sampleTime/2*1000, noise(:,1:N), 'k-');
    for k = 1:numel(hPlot)
       hPlot(k).Color(4) = 0.1;  % 5% transparent
    end
    % Last instance in green
    hPlot(numel(hPlot)).Color = [0 1 0 1];
    hPlot(numel(hPlot)).LineWidth = 1.5;
    
    set(gca, 'YLim', visualizeAmplitudeRange*[-1 1], 'XTick', [0:25:500]);
    xlabel('\it time (msec)');
    ylabel('\it amplitude (pAmps)');
    set(gca, 'FontSize', 14);
    axis 'square'; grid on; box on;
    % Add the noise amplitude histogram
    hold on;
    [NN, edges] = histcounts(noise, [-20:1:10]);
    barh(edges(1:(end-1)),NN/max(NN)*timeAxis(end)*0.5, 1, 'FaceAlpha', 0.4, 'FaceColor', [1 1 1]);

    % Export Fig
    NicePlot.exportFigToPDF('photocurrentNoiseModel.pdf', hFig, 300)
    
end


function [timeAxis, noise, noiseSpectrum] = generateNoise(nSamples, sampleTime, noiseSPD, nInstances)
    % Generate white Gaussian noise
    noise = randn(nInstances, nSamples);
    noiseFFT = fft(noise, [], 2);
    
    % Make-up the negative frequency part
    noiseSPDfull = [noiseSPD noiseSPD(end:-1:1)];
    noiseSPDfull = noiseSPDfull(1:nSamples);
    % Adjust the spectral power distribution of the noise
    noiseFFT = bsxfun(@times, noiseFFT, sqrt(noiseSPDfull));
    % Get the noise by ifft
    noise = real(ifft(noiseFFT, [], 2)) / sqrt(2 * sampleTime);
    timeAxis = (1:size(noise,2))*sampleTime;
    ft = fft(noise,[],2);
    noiseSpectrum = 2 * sampleTime/numel(timeAxis) * (ft .* conj(ft));
    % mean over instances
    noiseSpectrum = mean(noiseSpectrum, 1);
end

function [freqAxis, noiseSPD, noiseSPDlowFreq, ...
    noiseSPDhighFreq, cornerFreqLow, cornerFreqHigh] = ...
        generateSPD(nSamples, sampleTime)
    
    k = ceil((nSamples-1)/2);
    freqAxis = (0:k) / sampleTime / nSamples;
    LorentzCoeffs = [0.16  55  4;
                 0.045 190 2.5];
    noiseSPDlowFreq = lorentzSum(LorentzCoeffs(1,:), freqAxis);
    noiseSPDhighFreq = lorentzSum(LorentzCoeffs(2,:), freqAxis);
    noiseSPD = lorentzSum(LorentzCoeffs, freqAxis);
    
    [~,idxLow] = min(abs(noiseSPDlowFreq - noiseSPDlowFreq(1)/2));
    [~,idxHigh] = min(abs(noiseSPDhighFreq - noiseSPDhighFreq(1)/2));
    
    cornerFreqLow = freqAxis(idxLow);
    cornerFreqHigh = freqAxis(idxHigh);
end


