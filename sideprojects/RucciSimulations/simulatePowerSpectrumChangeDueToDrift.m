function simulatePowerSpectrumChangeDueToDrift

    % Generate eye movements
    [emPosArcMin, timeAxis] = generateFEM();
    plotFEM(emPosArcMin, timeAxis);
    
    % Make each pixel equal to 1 cone aperture
    coneApertureMicrons = 1.5;
    micronsPerDegree = 300;
    pixelSizeArcMin = 1*coneApertureMicrons / micronsPerDegree * 60;
    
    % Stimulus size
    stimSizeArcMin = 400;
    
    instancesNum = 25;
    
    % Grating params
    gratingParams.oriDegs = 45;
    gratingParams.sigmaArcMin = 12;
    gratingParams.contrast = 0.2;
    
    % Noise params
    % Noise params for low frequency stimulus
    gratingParams.sfCPD = 4;
    noiseParams.spectrumShape = 'highPassCornerFreq';
    noiseParams.cornerFrequencyCPD = 10;  % relevant for 'lowPassCornerFreq' and 'lhighPassCornerFreq'
    noiseParams.steepness = 2;
    % Generate stimulus
    [lowFrequencyStimulusImage, spatialSupportDegs, lowFrequencySpectrum, spatialFrequencySupport] = ...
        generateStimulusImage(stimSizeArcMin, pixelSizeArcMin,  gratingParams, noiseParams, instancesNum);
    
    % Plot stimulus
    plotStimulusAndSpectrum(lowFrequencyStimulusImage, spatialSupportDegs, lowFrequencySpectrum, spatialFrequencySupport, 1);
    
    % Noise params for highfrequency stimulus
    gratingParams.sfCPD = 11;
    noiseParams.spectrumShape = 'lowPassCornerFreq';
    noiseParams.cornerFrequencyCPD = 5;  % relevant for 'lowPassCornerFreq' and 'lhighPassCornerFreq'
    % Generate stimulus
    [highFrequencyStimulusImage, spatialSupportDegs, highFrequencySpectrum, spatialFrequencySupport] = ...
        generateStimulusImage(stimSizeArcMin, pixelSizeArcMin,  gratingParams, noiseParams, instancesNum);
    
    % Plot stimulus
    plotStimulusAndSpectrum(highFrequencyStimulusImage, spatialSupportDegs, highFrequencySpectrum, spatialFrequencySupport, 2);
    
    
    % Noise params for zero stimulus
    gratingParams.sfCPD = 0;
    noiseParams.spectrumShape = '1overF';
    % Generate stimulus
    [highFrequencyStimulusImage, spatialSupportDegs, highFrequencySpectrum, spatialFrequencySupport] = ...
        generateStimulusImage(stimSizeArcMin, pixelSizeArcMin,  gratingParams, noiseParams, instancesNum);
    
    % Plot stimulus
    plotStimulusAndSpectrum(highFrequencyStimulusImage, spatialSupportDegs, highFrequencySpectrum, spatialFrequencySupport, 3);
    
    
    % Noise params for zero stimulus
    gratingParams.sfCPD = 11;
    noiseParams.spectrumShape = 'none';
    % Generate stimulus
    [highFrequencyStimulusImage, spatialSupportDegs, highFrequencySpectrum, spatialFrequencySupport] = ...
        generateStimulusImage(stimSizeArcMin, pixelSizeArcMin,  gratingParams, noiseParams, instancesNum);
    
    % Plot stimulus
    plotStimulusAndSpectrum(highFrequencyStimulusImage, spatialSupportDegs, highFrequencySpectrum, spatialFrequencySupport, 4);
    
    
end

% Generate noise with spectrum f^{alpha}
function [stimulusImage, spatialSupportDegs, spectrum, spatialFrequencySupport] = ...
    generateStimulusImage(stimSizeArcMin, pixelSizeArcMin,  gratingParams, noiseParams, instancesNum)
    N = round(stimSizeArcMin / pixelSizeArcMin);
    if (mod(N,2) == 1)
        N = N + 1;
    end
    
    % Generate spatial support
    xAxis = (1:N)*pixelSizeArcMin;
    xAxis = xAxis - mean(xAxis);
    spatialSupportDegs = xAxis / 60;
    
    % Generate grating component
    [gratingComponent,gaussianEnvelope] = generateGratingComponent(spatialSupportDegs, gratingParams);
    
    % Generate noise component
    for instanceNo = 1:instancesNum
        noiseComponent = generateNoiseComponent(spatialSupportDegs, noiseParams);
 
        % Add components
        stimulus = gaussianEnvelope .* (noiseComponent + gratingComponent);
    
        % Normalize to [-1 1]
        stimulusNorm = stimulus / max(abs(stimulus(:)));
    
        % Scale to [0 .. 1]
        stimulusImage(instanceNo,:,:) = 0.5*(1+stimulusNorm);
    
        % Compute stimulus spectrum and spatialFrequencySupport
        spectrum(instanceNo,:,:) = fftshift(abs(fft2(squeeze(stimulusImage(instanceNo,:,:)))));
        positiveSFindices = (N/2+1):N;
        spacingDegs = pixelSizeArcMin/60;
        sfMaxCPD = 1/(2*spacingDegs);
        spatialFrequencySupport = (0:(numel(positiveSFindices)-1))/(numel(positiveSFindices)-1) * sfMaxCPD;
    end
    
end


function noiseImage = generateNoiseComponent(spatialSupportDegs, noiseParams)
    N = length(spatialSupportDegs);
    
    % Determine max SF
    spacingDegs = spatialSupportDegs(2)-spatialSupportDegs(1);
    sfMaxCPD = 1/(2*spacingDegs);
    
    % Generate the grid of frequencies. First quadrant are positive frequencies.  
    % Zero frequency is at u(1,1).
    sf = [(0:floor(N/2)) -(ceil(N/2)-1:-1:1)]'/N;
    % Replicate frequencies
    sfX = repmat(sf,1,N);

    % Generate the normalized spatial frequency grid (fftshifted)
    normalizedSpatialFrequencyGridFFTshift = sqrt(sfX.^2 + (sfX').^2);
    
    
    % Generate oneOverFspectrum 1/(f^alpha)
    oneOverFspectrum = normalizedSpatialFrequencyGridFFTshift .^(-1.0);
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
            xTerm = normalizedSpatialFrequencyGridFFTshift.^noiseParams.steepness;
            x50Term = (noiseParams.cornerFrequencyCPD/sfMaxCPD)^noiseParams.steepness;
            noiseSpectrumFilter = xTerm ./ (xTerm + x50Term);
             
            if (strcmp(noiseParams.spectrumShape,'lowPassCornerFreq'))
                noiseSpectrumFilter = 1-noiseSpectrumFilter;
            end
            
            % Apply filter to 1/F noise
            noiseSpectrum = oneOverFspectrum .* noiseSpectrumFilter;
    end
    

    % Generate a grid of random phase shifts
    phi = rand(N);

    % Inverse Fourier transform to obtain the the spatial pattern
    noiseImage = ifft2(noiseSpectrum .* (cos(2*pi*phi)+1i*sin(2*pi*phi)));

    % Pick the real component
    noiseImage = real(noiseImage);
    noiseImage = noiseImage / max(abs(noiseImage(:)));
    
    if (strcmp(noiseParams.spectrumShape,'none'))
        disp('all zeros')
        noiseImage = zeros(size(noiseImage));
    end
     
    plotNoiseSpectrum = false;
    if (plotNoiseSpectrum)
        figure(222);
        clf;
        subplot(2,1,1);
        sfXX = fftshift(sfX);
        spatialFrequencySupportCPD = sfXX(2:end-1,1) * sfMaxCPD;
        imagesc(spatialFrequencySupportCPD, spatialFrequencySupportCPD, log10(fftshift(noiseSpectrum/max(abs(noiseSpectrum(:))))));
        axis 'image'
        set(gca, 'XLim', 20*[-1 1], 'YLim',20*[-1 1]);

        subplot(2,1,2);
        spatialFrequencySupportCPD = sfX(1:N/2,1) * sfMaxCPD;
        plot(spatialFrequencySupportCPD, log10(noiseSpectrum(1:N/2,1)), 'k-');
        set(gca, 'XLim', 20*[-1 1]);
        colormap(gray);
        axis 'square';
    end
    
end

function [gratingImage, gaussianEnvelope] = generateGratingComponent(spatialSupportDegs, gratingParams)

    N = length(spatialSupportDegs);
    [xx,yy] = meshgrid(spatialSupportDegs, spatialSupportDegs);
    gratingParams.sfCPD
    gratingImage = gratingParams.contrast * cos(2*pi*gratingParams.sfCPD*(xx*cosd(gratingParams.oriDegs) + yy*sind(gratingParams.oriDegs)));
    gaussianEnvelope = exp(-0.5*((xx/(gratingParams.sigmaArcMin/60)).^2)) .* exp(-0.5*((yy/(gratingParams.sigmaArcMin/60)).^2));
end


function plotStimulusAndSpectrum(stimulusImage, spatialSupportDegs, spectrum, spatialFrequencySupport, figNo)
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [200 100+(figNo-1)*600 1500 600], 'Color', [1 1 1]);
    colormap(gray(1024));
    
    instancesNum = size(spectrum,1);
    N = size(spectrum,2);
    positiveSFindices = (N/2+1):N;
    diagonalFrequencySupport = sqrt(spatialFrequencySupport.^2+spatialFrequencySupport.^2);
    spatialSupportArcMin = spatialSupportDegs*60;
    
    visualizedInstancesNum = 5;
    for instanceNo = 1:visualizedInstancesNum
        subplot(2,visualizedInstancesNum+1,instanceNo);
        
        imagesc(spatialSupportArcMin, spatialSupportArcMin,  (squeeze(stimulusImage(instanceNo,:,:))).^0.5);
        axis 'image';
        set(gca, 'CLim', [0 1], 'FontSize', 14, 'XLim', 60*[-1 1], 'YLim', 60*[-1 1]);
        xlabel('space (arc min)');
        if (instanceNo == 1)
            ylabel('space (arc min)');
        end
        title(sprintf('instance %2.0f', instanceNo));
        
        subplot(2,visualizedInstancesNum+1,instanceNo+visualizedInstancesNum+1);
        % 0 deg slice
        spectrumSlice1 = squeeze(spectrum(instanceNo,N/2+1,positiveSFindices));
        % 45 deg slice
        diagonalSlice = diag(squeeze(spectrum(instanceNo,:,:)));
        spectrumSlice2 = squeeze(diagonalSlice(positiveSFindices));
        % 90 deg slice
        spectrumSlice3 = squeeze(spectrum(instanceNo,positiveSFindices,N/2+1));
        plot(spatialFrequencySupport, 10*log10(spectrumSlice1), 'k-', 'LineWidth', 1.5); hold on
        plot(diagonalFrequencySupport, 10*log10(spectrumSlice2), 'b-', 'LineWidth', 1.5);
        plot(spatialFrequencySupport, 10*log10(spectrumSlice3), 'r-', 'LineWidth', 1.5);
        hl = legend({'0', '45', '90'});
        set(gca, 'XScale', 'log', 'XLim', [0.5 60], 'YLim', [0 40]);
        xlabel('spacial frequency (c/deg)');
        if (instanceNo == 1)
            ylabel('power (dB)');
        end
        grid on
        set(gca, 'XTick', [0.5 1 2 3 5  10 20 50  100], 'FontSize', 14);
        drawnow
    end
    
    spectrum = squeeze(mean(spectrum,1));
    % 0 deg slice
    spectrumSlice1Mean = squeeze(spectrum(N/2+1,positiveSFindices));
    % 45 deg slice
    diagonalSlice = diag(squeeze(spectrum(:,:)));
    spectrumSlice2Mean = squeeze(diagonalSlice(positiveSFindices));
    % 90 deg slice
    spectrumSlice3Mean = squeeze(spectrum(positiveSFindices,N/2+1));
        
    subplot(2,visualizedInstancesNum+1,visualizedInstancesNum*2+2);
    % Mean spectra across instances
    plot(spatialFrequencySupport, 10*log10(spectrumSlice1Mean), 'k-', 'LineWidth', 1.5); hold on
    plot(diagonalFrequencySupport, 10*log10(spectrumSlice2Mean), 'b-', 'LineWidth', 1.5);
    plot(spatialFrequencySupport, 10*log10(spectrumSlice3Mean), 'r-', 'LineWidth', 1.5);
    hl = legend({'0', '45', '90'});
    set(gca, 'XScale', 'log', 'XLim', [0.5 60], 'YLim', [0 40]);
    xlabel('spacial frequency (c/deg)');
    if (instanceNo == 1)
        ylabel('power (dB)');
    end
    grid on
    title('mean spectra');
    set(gca, 'XTick', [0.5 1 2 3 5 10 20 50  100], 'FontSize', 14);
    drawnow
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
    
    hFig = figure(1000); clf;
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
