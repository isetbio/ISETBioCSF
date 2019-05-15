function simulatePowerSpectrumChangeDueToDrift

    % Make each pixel equal to 1 cone aperture
    coneApertureMicrons = 1.5;
    micronsPerDegree = 300;
    pixelSizeArcMin = coneApertureMicrons / micronsPerDegree * 60;
    
    % Generate stimulus
    stimRadiusArcMin = 30;
    stimSizeArcMin = 100;
    
    [stimulusImage, spatialSupportDegs, spectrum, spatialFrequencySupport] = ...
        generateStimulusImage(stimSizeArcMin, pixelSizeArcMin, stimRadiusArcMin, -1);
    plotStimulusAndSpectrum(stimulusImage, spatialSupportDegs, spectrum, spatialFrequencySupport);
    
    % Generate eye movements
    [emPosArcMin, timeAxis] = generateFEM();
    plotFEM(emPosArcMin, timeAxis);
    
end

% Generate noise with spectrum f^{alpha}
function [stimulusImage, spatialSupportDegs, spectrum, spatialFrequencySupport] = generateStimulusImage(stimSizeArcMin, pixelSizeArcMin, stimRadiusArcMin, alpha)
    N = round(stimSizeArcMin / pixelSizeArcMin);
    if (mod(N,2) == 1)
        N = N + 1;
    end
    
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
    x = ifft2(noiseSpectrum .* (cos(2*pi*phi)+1i*sin(2*pi*phi)));

    % Pick the real component
    x = real(x);
    
    % Scale to [0 .. 1]
    stimulusImage = 0.5*(1+x / max(abs(x(:))));
    
    % Generate spatial support
    xAxis = (1:N)*pixelSizeArcMin;
    xAxis = xAxis - mean(xAxis);
    spatialSupportDegs = xAxis / 60;

    % Mask stimulus to 0.5 outside of stimRadiusArcMin
    [xx,yy] = meshgrid(spatialSupportDegs, spatialSupportDegs);
    r = sqrt(xx.^2 + yy.^2);
    idx = find(r > stimRadiusArcMin/60);
    stimulusImage(idx) = 0.5;
    
    % Compute stimulus spectrum and spatialFrequencySupport
    spectrum = fftshift(abs(fft2(stimulusImage)));
    normalizedSpatialFrequencyGrid = ifftshift(normalizedSpatialFrequencyGridFFTshift);
    positiveSFindices = (N/2+1):N;
    spacingDegs = pixelSizeArcMin/60;
    sfMaxCPD = 1/(2*spacingDegs);
    spatialFrequencySupport = normalizedSpatialFrequencyGrid(positiveSFindices, N/2+1) * sfMaxCPD;
end

function plotStimulusAndSpectrum(stimulusImage, spatialSupportDegs, spectrum, spatialFrequencySupport)
    figure(1); clf;
    subplot(1,2,1);
    spatialSupportArcMin = spatialSupportDegs*60;
    imagesc(spatialSupportArcMin, spatialSupportArcMin,  stimulusImage.^0.5);
    set(gca, 'CLim', [0 1]);
    xlabel('space (arc min)');
    axis 'image';
    colormap(gray);

    subplot(1,2,2);
    N = size(spectrum,1);
    positiveSFindices = (N/2+1):N;
    spectrumSlice1 = squeeze(spectrum(positiveSFindices,N/2+1));
    spectrumSlice2 = squeeze(spectrum(N/2+1,positiveSFindices));
    meanSpectrumSlice = 0.5*(spectrumSlice1(:)+spectrumSlice2(:));
    plot(spatialFrequencySupport, 10*log10(meanSpectrumSlice), 'k-', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log');
    axis 'square';
    xlabel('spacial frequency (c/deg)');
    ylabel('power (dB)');
    set(gca, 'XTick', [0.5 1 2 4 8 16 32 64]);
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
    
    figure(2); clf;
    subplot(1,2,1);
    plot(timeAxis, xPosArcMin, 'r-', 'LineWidth', 1.5); hold on
    plot(timeAxis, yPosArcMin, 'b-', 'LineWidth', 1.5);
    set(gca, 'YLim', posRange);
    xlabel('time (msec)');
    ylabel('space (arc min)');
    axis 'square';
    
    subplot(1,2,2);
    plot(xPosArcMin, yPosArcMin, 'k-', 'LineWidth', 1.5);
    set(gca, 'XLim', posRange, 'YLim', posRange);
    axis 'square';
    xlabel('space (arc min)');
end
