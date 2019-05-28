% Generate noise with spectrum f^{alpha}
function [stimulus, noiseNorm] = generateStimulusImage(stimSizeArcMin, pixelSizeArcMin,  gratingParams, noiseParams, noiseNorm, instancesNum, figNo)

    N = round(stimSizeArcMin / pixelSizeArcMin);
    if (mod(N,2) == 1)
        N = N + 1;
    end
    
    % Generate spatial support
    xAxis = (1:N)*pixelSizeArcMin;
    xAxis = xAxis - mean(xAxis);
    spatialSupportDegs = xAxis / 60;
    
    % Generate spatial frequency support
    positiveSFindices = (N/2+1):N;
    spacingDegs = pixelSizeArcMin/60;
    sfMaxCPD = 1/(2*spacingDegs);
    spatialFrequencySupport = (0:(numel(positiveSFindices)-1))/(numel(positiveSFindices)-1) * sfMaxCPD;
        
    % Generate grating component
    [gratingComponent,gaussianEnvelope] = generateGratingComponent(spatialSupportDegs, gratingParams);
    
    computePower = false;
    if (computePower)
        % Compute power
        gratingImagePower = computeImageTotalPower(gratingComponent);
    end
    
    % Generate as many noise components as there are instances
    for instanceNo = 1:instancesNum
        [noiseComponent, noiseNorm] = generateNoiseComponent(spatialSupportDegs, noiseParams, noiseNorm);
 
        % Compute total noise power
        if (computePower)
            noiseImagePower(instanceNo) = computeImageTotalPower(noiseComponent);
        end
        
        % Add components
        stimulus = gaussianEnvelope .* (noiseComponent + gratingComponent);
    
        % Normalize to [-1 1]
        stimulusImage(instanceNo,:,:) = stimulus / max(abs(stimulus(:)));
    
        % Compute stimulus spectrum and spatialFrequencySupport
        spectrum(instanceNo,:,:) = fftshift(abs(fft2(squeeze(stimulusImage(instanceNo,:,:)))));
    end
    
    if (computePower)
        fprintf('Noise image power: %f\n', mean(noiseImagePower));
        fprintf('grating image power: %f\n', gratingImagePower);
    end
    
    stimulus = struct(...
        'image', stimulusImage, ...
        'spatialSupportDegs', spatialSupportDegs, ...
        'spectrum', spectrum, ...
        'spatialFrequencySupport', spatialFrequencySupport);
    
    plotStimulusAndSpectrum(stimulus, figNo);
    
end

function [noiseImage, noiseNorm] = generateNoiseComponent(spatialSupportDegs, noiseParams, noiseNorm)
    N = length(spatialSupportDegs);
    
    % Determine max SF
    spacingDegs = spatialSupportDegs(2)-spatialSupportDegs(1);
    sfMaxCPD = 1/(spacingDegs);
    
    % Generate the grid of frequencies. First quadrant are positive frequencies.  
    % Zero frequency is at u(1,1).
    sf = [(0:floor(N/2)) -(ceil(N/2)-1:-1:1)]'/N;
    % Replicate frequencies
    sfX = repmat(sf,1,N);

    % Generate the normalized spatial frequency grid (fftshifted)
    normalizedSpatialFrequencyGridFFTshift = sqrt(sfX.^2 + (sfX').^2);
    
    
    % Generate oneOverFspectrum 1/(f^alpha)
    alpha = 1;
    oneOverFspectrum = 1 ./ (normalizedSpatialFrequencyGridFFTshift.^alpha);
    % Set any infinities to zero
    oneOverFspectrum(oneOverFspectrum==inf) = 0;
    oneOverFspectrum = oneOverFspectrum / max(oneOverFspectrum(:));
    oneOverFspectrum(1,1) = 1;
            
    % Generate the spectrum shape
    switch(noiseParams.spectrumShape)
        case '1overF'
            noiseSpectrum = oneOverFspectrum;
            
        case 'none'
            noiseSpectrum = oneOverFspectrum*0;
            
        case {'lowPassCornerFreq','highPassCornerFreq'}
            
            % generate low-pass or high-pass filter
            a = (normalizedSpatialFrequencyGridFFTshift).^noiseParams.steepness;
            b = (noiseParams.cornerFrequencyCPD/sfMaxCPD)^noiseParams.steepness;
            frequencyFilter = a ./ (a + b);
            
            if (strcmp(noiseParams.spectrumShape,'lowPassCornerFreq'))
                frequencyFilter = 1-frequencyFilter;
            end
            
            % Apply filter to 1/F noise
            noiseSpectrum = oneOverFspectrum .* frequencyFilter;
    end
    
    % Generate a grid of random phases
    phaseComponent = 2*pi*rand(size(noiseSpectrum));

    % Generate full FFT
    noiseFFT = noiseSpectrum .* (cos(phaseComponent)+1i*sin(phaseComponent));
    
    % Inverse Fourier transform to obtain the the spatial pattern
    noiseImage = ifft2(noiseFFT);

    % Pick the real component
    noiseImage = real(noiseImage);
    if (isnan(noiseNorm))
        noiseNorm = sum(sum(noiseImage.^2));
    end
    noiseImage = noiseImage / noiseNorm;
    
    if (strcmp(noiseParams.spectrumShape,'none'))
        noiseImage = zeros(size(noiseImage));
    end
     
    % Plot the noise spectrum for debug purposes
    plotNoiseSpectrum = false;
    if (plotNoiseSpectrum)
        figure(223); clf;
        maxSFdisplayed = 20;
        sfXX = fftshift(sfX);
        spatialFrequencySupportCPD = sfXX(2:end,1) * sfMaxCPD;

        subplot(2,3,1);
        m = fftshift(oneOverFspectrum);
        imagesc(spatialFrequencySupportCPD, spatialFrequencySupportCPD, m);
        axis 'image';
        set(gca, 'XLim', maxSFdisplayed*[-1 1], 'YLim', maxSFdisplayed*[-1 1]);
        title('1/F spectrum');
        subplot(2,3,4);
        plot(spatialFrequencySupportCPD, m(size(m,1)/2+1,2:end), 'rs-')
        set(gca, 'YScale', 'log', 'YLim', 10.^[-3 0], 'XLim', maxSFdisplayed*[-1 1]);

        subplot(2,3,2);
        m = fftshift(frequencyFilter);
        imagesc(spatialFrequencySupportCPD, spatialFrequencySupportCPD, m);
        axis 'image';
        set(gca, 'XLim', maxSFdisplayed*[-1 1], 'YLim', maxSFdisplayed*[-1 1]);

        title('frequency filter');
        subplot(2,3,5);
        plot(spatialFrequencySupportCPD, m(size(m,1)/2+1,2:end), 'rs-')
        set(gca, 'YScale', 'log', 'YLim', 10.^[-3 0], 'XLim', maxSFdisplayed*[-1 1]);

        subplot(2,3,3);
        m = fftshift(noiseSpectrum);
        imagesc(spatialFrequencySupportCPD, spatialFrequencySupportCPD, m);
        axis 'image';
        set(gca, 'XLim', maxSFdisplayed*[-1 1], 'YLim', maxSFdisplayed*[-1 1]);

        title('noise spectrum');
        subplot(2,3,6);
        plot(spatialFrequencySupportCPD, m(size(m,1)/2+1,2:end), 'rs-');
        hold on;
        m2 = abs(fftshift(fft2(noiseImage)));
        m2 = m2 / max(m2(:)) * max(m(:));
        plot(spatialFrequencySupportCPD, m2(size(m,1)/2+1,2:end), 'bs-');
        set(gca, 'YScale', 'log', 'YLim', 10.^[-3 0], 'XLim', maxSFdisplayed*[-1 1]);
        colormap(gray);
        pause
    end

end

function noiseImagePowerDB = computeImageTotalPower(img)
     noiseImagePowerDB = 10*log10(sum(sum((abs(fft2(img))).^2)));
end

function [gratingImage, gaussianEnvelope] = generateGratingComponent(spatialSupportDegs, gratingParams)
    N = length(spatialSupportDegs);
    [xx,yy] = meshgrid(spatialSupportDegs, spatialSupportDegs);
    gratingImage = gratingParams.contrast * cos(2*pi*gratingParams.sfCPD*(xx*cosd(gratingParams.oriDegs) + yy*sind(gratingParams.oriDegs)));
    gaussianEnvelope = exp(-0.5*((xx/(gratingParams.sigmaArcMin/60)).^2)) .* exp(-0.5*((yy/(gratingParams.sigmaArcMin/60)).^2));
end

function plotStimulusAndSpectrum(stimulus, figNo)
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [200 100+(figNo-1)*600 1500 600], 'Color', [1 1 1]);
    colormap(gray(1024));
    
    N = size(stimulus.spectrum,2);
    positiveSFindices = (N/2+1):N;
    diagonalFrequencySupport = sqrt(stimulus.spatialFrequencySupport.^2+stimulus.spatialFrequencySupport.^2);
    spatialSupportArcMin = stimulus.spatialSupportDegs*60;
    
    visualizedInstancesNum = min([5 size(stimulus.image,1)]);
    xRange = max(spatialSupportArcMin(:))*[-1 1];
    for instanceNo = 1:visualizedInstancesNum
        
        % Plot the stimulus
        subplot(2,visualizedInstancesNum+1,instanceNo);
        % exponent of 0.5 to gamma correct
        imagesc(spatialSupportArcMin, spatialSupportArcMin,  (0.5*(1+squeeze(stimulus.image(instanceNo,:,:)))).^0.5);
        
        axis 'image';
        set(gca, 'CLim', [0 1], 'FontSize', 14, 'XLim', xRange, 'YLim', xRange);
        xlabel('space (arc min)');
        if (instanceNo == 1)
            ylabel('space (arc min)');
        end
        title(sprintf('instance %2.0f', instanceNo));
        
        % Plot slices (different orientations) through the 2D stimulus spectrum
        subplot(2,visualizedInstancesNum+1,instanceNo+visualizedInstancesNum+1);
        % 0 deg slice
        spectrumSlice1 = squeeze(stimulus.spectrum(instanceNo,N/2+1,positiveSFindices));
        % 45 deg slice
        diagonalSlice = diag(squeeze(stimulus.spectrum(instanceNo,:,:)));
        spectrumSlice2 = squeeze(diagonalSlice(positiveSFindices));
        % 90 deg slice
        spectrumSlice3 = squeeze(stimulus.spectrum(instanceNo,positiveSFindices,N/2+1));
        plot(stimulus.spatialFrequencySupport, 10*log10(spectrumSlice1), 'k-', 'LineWidth', 1.5); hold on
        plot(diagonalFrequencySupport, 10*log10(spectrumSlice2), 'b-', 'LineWidth', 1.5);
        plot(stimulus.spatialFrequencySupport, 10*log10(spectrumSlice3), 'r-', 'LineWidth', 1.5);
        hl = legend({'0', '45', '90'});
        set(gca, 'XScale', 'log', 'XLim', [0.5 100], 'YLim', [0 40]);
        xlabel('spatial frequency (c/deg)');
        if (instanceNo == 1)
            ylabel('power (dB)');
        end
        grid on
        set(gca, 'XTick', [0.5 1 2 3 5  10 20 50  100], 'FontSize', 14);
        drawnow
    end
    
    % Compute the mean spectrum over all instances
    stimulus.spectrum = squeeze(mean(stimulus.spectrum,1));
    % 0 deg slice through the mean spectrum
    spectrumSlice1Mean = squeeze(stimulus.spectrum(N/2+1,positiveSFindices));
    % 45 deg slice through the mean spectrum
    diagonalSlice = diag(squeeze(stimulus.spectrum(:,:)));
    spectrumSlice2Mean = squeeze(diagonalSlice(positiveSFindices));
    % 90 deg slice through the mean spectrum
    spectrumSlice3Mean = squeeze(stimulus.spectrum(positiveSFindices,N/2+1));
        
    subplot(2,visualizedInstancesNum+1,visualizedInstancesNum*2+2);
    % Mean spectra across instances
    plot(stimulus.spatialFrequencySupport, 10*log10(spectrumSlice1Mean), 'k-', 'LineWidth', 1.5); hold on
    plot(diagonalFrequencySupport, 10*log10(spectrumSlice2Mean), 'b-', 'LineWidth', 1.5);
    plot(stimulus.spatialFrequencySupport, 10*log10(spectrumSlice3Mean), 'r-', 'LineWidth', 1.5);
    hl = legend({'0', '45', '90'});
    set(gca, 'XScale', 'log', 'XLim', [0.5 60], 'YLim', [0 40]);
    xlabel('spatial frequency (c/deg)');
    if (instanceNo == 1)
        ylabel('power (dB)');
    end
    grid on
    title('mean spectra');
    set(gca, 'XTick', [0.5 1 2 3 5 10 20 50  100], 'FontSize', 14);
    drawnow
end

