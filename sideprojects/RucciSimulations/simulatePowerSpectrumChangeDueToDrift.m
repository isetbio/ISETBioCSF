function simulatePowerSpectrumChangeDueToDrift

    % Generate eye movements
    [emPosArcMin, timeAxis] = generateFEM();
    plotFEM(emPosArcMin, timeAxis);
    
    % Make each pixel equal to 1 cone aperture
    coneApertureMicrons = 1.5;
    micronsPerDegree = 300;
    pixelSizeArcMin = coneApertureMicrons / micronsPerDegree * 60;
    
    % Stimulus size
    stimSizeArcMin = 120;
    
    % Stimulus radius
    stimRadiusArcMin = 40;
        
    % Grating params
    gratingParams.sfCPD = 15;
    gratingParams.oriDegs = 85;
    gratingParams.sigmaArcMin = 10;
    gratingParams.contrast = 0.7;
    
    % Noise params
    alpha = -1;  % noise spectrum: sf^-alpha
    
    % Generate stimulus
    [stimulusImage, spatialSupportDegs, spectrum, spatialFrequencySupport] = ...
        generateStimulusImage(stimSizeArcMin, pixelSizeArcMin, stimRadiusArcMin, gratingParams, alpha);
    
    % Plot stimulus
    plotStimulusAndSpectrum(stimulusImage, spatialSupportDegs, spectrum, spatialFrequencySupport);
end

% Generate noise with spectrum f^{alpha}
function [stimulusImage, spatialSupportDegs, spectrum, spatialFrequencySupport] = generateStimulusImage(stimSizeArcMin, pixelSizeArcMin, stimRadiusArcMin, gratingParams, alpha)
    N = round(stimSizeArcMin / pixelSizeArcMin);
    if (mod(N,2) == 1)
        N = N + 1;
    end
    
    % Generate spatial support
    xAxis = (1:N)*pixelSizeArcMin;
    xAxis = xAxis - mean(xAxis);
    spatialSupportDegs = xAxis / 60;
    [xx,yy] = meshgrid(spatialSupportDegs, spatialSupportDegs);
    
    % Generate the grid of frequencies. First quadrant are positive frequencies.  
    % Zero frequency is at u(1,1).
    sf = [(0:floor(N/2)) -(ceil(N/2)-1:-1:1)]'/N;
    % Reproduce these frequencies
    sfX = repmat(sf,1,N);

    % Generate the power spectrum
    normalizedSpatialFrequencyGridFFTshift = sqrt(sfX.^2 + (sfX').^2);
    noiseSpectrum = normalizedSpatialFrequencyGridFFTshift.^alpha;

    % Set any infinities to zero
    noiseSpectrum(noiseSpectrum==inf) = 0;

    % Generate a grid of random phase shifts
    phi = rand(N);

    % Inverse Fourier transform to obtain the the spatial pattern
    noiseImage = ifft2(noiseSpectrum .* (cos(2*pi*phi)+1i*sin(2*pi*phi)));

    % Pick the real component
    noiseImage = real(noiseImage);
    noiseImage = noiseImage / max(abs(noiseImage(:)));
    
    % Add the grating
    gratingImage = gratingParams.contrast * sin(2*pi*gratingParams.sfCPD*(xx*cosd(gratingParams.oriDegs) + yy*sind(gratingParams.oriDegs)));
    gaussianEnvelope = exp(-0.5*((xx/(gratingParams.sigmaArcMin/60)).^2)) .* exp(-0.5*((yy/(gratingParams.sigmaArcMin/60)).^2));
    gaborImage = gratingImage .* gaussianEnvelope;
    noiseImage = 0.5*(noiseImage + gaborImage);
    
    % Scale to [0 .. 1]
    stimulusImage = 0.5*(1+noiseImage / max(abs(noiseImage(:))));
    
    % Mask stimulus to 0.5 outside of stimRadiusArcMin
    r = sqrt(xx.^2 + yy.^2);
    idx = find(r > stimRadiusArcMin/60);
    stimulusImage(idx) = 0.5;
    
    % Compute stimulus spectrum and spatialFrequencySupport
    spectrum = fftshift(abs(fft2(stimulusImage)));
    normalizedSpatialFrequencyGrid = ifftshift(normalizedSpatialFrequencyGridFFTshift);
    positiveSFindices = (N/2+1):N;
    spacingDegs = pixelSizeArcMin/60;
    sfMaxCPD = 1/spacingDegs;
    spatialFrequencySupport = normalizedSpatialFrequencyGrid(positiveSFindices, N/2+1) * sfMaxCPD;
end

function plotStimulusAndSpectrum(stimulusImage, spatialSupportDegs, spectrum, spatialFrequencySupport)
    hFig = figure(1); clf;
    set(hFig, 'Position', [500 500 900 400], 'Color', [1 1 1]);
    
    subplot(1,2,1);
    spatialSupportArcMin = spatialSupportDegs*60;
    imagesc(spatialSupportArcMin, spatialSupportArcMin,  stimulusImage.^0.5);
    set(gca, 'CLim', [0 1], 'FontSize', 14);
    xlabel('space (arc min)');
    ylabel('space (arc min)');
    axis 'image';
    colormap(gray);

    subplot(1,2,2);
    N = size(spectrum,1);
    positiveSFindices = (N/2+1):N;
    % 0 deg slice
    spectrumSlice1 = squeeze(spectrum(N/2+1,positiveSFindices));
    % 45 deg slice
    diagonalSlice = diag(spectrum);
    diagonalFrequencySupport = sqrt(spatialFrequencySupport.^2+spatialFrequencySupport.^2);
    spectrumSlice2 = squeeze(diagonalSlice(positiveSFindices));
    % 90 deg slice
    spectrumSlice3 = squeeze(spectrum(positiveSFindices,N/2+1));
    plot(spatialFrequencySupport, 10*log10(spectrumSlice1), 'r-', 'LineWidth', 1.5); hold on
    plot(diagonalFrequencySupport, 10*log10(spectrumSlice2), 'b-', 'LineWidth', 1.5);
    plot(spatialFrequencySupport, 10*log10(spectrumSlice3), 'm-', 'LineWidth', 1.5);
    hl = legend({'0', '45', '90'});
    set(gca, 'XScale', 'log', 'XLim', [0.5 100], 'YLim', [0 40]);
    axis 'square';
    xlabel('spacial frequency (c/deg)');
    ylabel('power (dB)');
    grid on
    set(gca, 'XTick', [0.7 1 2 3 5 7 10 20 30 50 70 100], 'FontSize', 14);
end

function [emPosArcMin, timeAxis] = generateFEM()
    emDurationSeconds = 1;
    sampleTimeSeconds = 1.0 / 1000;
    nTrials = 64;
    
    microSaccadeType = 'none'; % , 'heatmap/fixation based', 'none'
    fixEMobj = fixationalEM();
    fixEMobj.randomSeed = 1;
    fixEMobj.microSaccadeType = microSaccadeType;
    
    computeVelocity = ~true;
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, ...
            nTrials, computeVelocity, 'useParfor', true);
        
    timeAxis = fixEMobj.timeAxis;
    emPosArcMin = fixEMobj.emPosArcMin;
    
end

function plotFEM(emPosArcMin, timeAxis)
    trialNo = 1;
    xPosArcMin = squeeze(emPosArcMin(trialNo,:,1));
    yPosArcMin = squeeze(emPosArcMin(trialNo,:,2));
    xPosArcMin = xPosArcMin - mean(xPosArcMin);
    yPosArcMin = yPosArcMin - mean(yPosArcMin);
    posRange = max(abs(emPosArcMin(:)))*[-1 1];
    
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 900 400], 'Color', [1 1 1]);
    subplot(1,2,1);
    plot(timeAxis*1000, xPosArcMin, 'r-', 'LineWidth', 1.5); hold on
    plot(timeAxis*1000, yPosArcMin, 'b-', 'LineWidth', 1.5);
    set(gca, 'YLim', posRange, 'FontSize', 14, 'XTick', 0:200:1000, 'YTick', -20:5:20);
    grid on
    xlabel('time (msec)');
    ylabel('space (arc min)');
    axis 'square';
    
    subplot(1,2,2);
    plot(xPosArcMin, yPosArcMin, 'k-', 'LineWidth', 1.5);
    set(gca, 'XLim', posRange, 'YLim', posRange, 'FontSize', 14, 'XTick', -20:5:20, 'YTick', -20:5:20);
    grid on
    axis 'square';
    xlabel('space (arc min)');
end
