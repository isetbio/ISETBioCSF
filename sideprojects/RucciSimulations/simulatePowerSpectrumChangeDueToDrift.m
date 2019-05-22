function simulatePowerSpectrumChangeDueToDrift

    % Generate eye movements
    nTrials = 10;
    emDurationSeconds = 0.2;
    [emPosArcMin, timeAxis] = generateFEM(nTrials, emDurationSeconds);
    
    
    % Plot the fixationalEM
    plotFEM(emPosArcMin, timeAxis);

    
    % Make each pixel equal to 1 cone aperture
    coneApertureMicrons = 2;
    micronsPerDegree = 300;
    pixelSizeArcMin = 1*coneApertureMicrons / micronsPerDegree * 60;
    
    % Stimulus size
    stimSizeArcMin = 150;
    
    % Grating params
    gratingParams.oriDegs = 90;
    gratingParams.sigmaArcMin = 12;
    gratingParams.contrast = 1;
    
    % Noise params
    % Noise params for low frequency stimulus
    gratingParams.sfCPD = 4;
    noiseParams.spectrumShape = 'highPassCornerFreq';
    noiseParams.cornerFrequencyCPD = 10;  % relevant for 'lowPassCornerFreq' and 'lhighPassCornerFreq'
    noiseParams.steepness = 60;
    % Generate stimulus
    [lowFrequencyStimulusImage, spatialSupportDegs, lowFrequencySpectrum, spatialFrequencySupport] = ...
        generateStimulusImage(stimSizeArcMin, pixelSizeArcMin,  gratingParams, noiseParams, nTrials);
    
    % Plot stimulus
    plotStimulusAndSpectrum(lowFrequencyStimulusImage, spatialSupportDegs, lowFrequencySpectrum, spatialFrequencySupport, 1);
    
    % Noise params for highfrequency stimulus
    gratingParams.sfCPD = 11;
    noiseParams.spectrumShape = 'lowPassCornerFreq';
    noiseParams.cornerFrequencyCPD = 5;  % relevant for 'lowPassCornerFreq' and 'lhighPassCornerFreq'
    % Generate stimulus
    [highFrequencyStimulusImage, spatialSupportDegs, highFrequencySpectrum, spatialFrequencySupport] = ...
        generateStimulusImage(stimSizeArcMin, pixelSizeArcMin,  gratingParams, noiseParams, nTrials);
    
    % Plot stimulus
    plotStimulusAndSpectrum(highFrequencyStimulusImage, spatialSupportDegs, highFrequencySpectrum, spatialFrequencySupport, 2);

    

    % Noise params for noise only stimulus
    gratingParams.contrast = 0.0;
    noiseParams.spectrumShape = '1overF';
    % Generate stimulus
    [noiseOnlyStimulusImage, spatialSupportDegs, noiseOnlySpectrum, spatialFrequencySupport] = ...
        generateStimulusImage(stimSizeArcMin, pixelSizeArcMin,  gratingParams, noiseParams, nTrials);
    
    % Plot stimulus
    plotStimulusAndSpectrum(noiseOnlyStimulusImage, spatialSupportDegs, noiseOnlySpectrum, spatialFrequencySupport, 3);
    
    
    % Noise params for grating only stimulus
    gratingParams.contrast = 1.0;
    noiseParams.spectrumShape = 'none';
    % Generate stimulus
    [gratingOnlyStimulusImage, spatialSupportDegs, gratingOnlySpectrum, spatialFrequencySupport] = ...
        generateStimulusImage(stimSizeArcMin, pixelSizeArcMin,  gratingParams, noiseParams, nTrials);
    
    % Plot stimulus
    plotStimulusAndSpectrum(gratingOnlyStimulusImage, spatialSupportDegs, gratingOnlySpectrum, spatialFrequencySupport, 4);
    
    
    %gratingStimulusSequence = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, gratingOnlyStimulusImage, spatialSupportDegs);
    highFrequencyStimulusSequence = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, highFrequencyStimulusImage, spatialSupportDegs);
    %lowFrequencyStimulusSequence = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, lowFrequencyStimulusImage, spatialSupportDegs);
    
end

function stimulusSequence = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, stimulusImage, spatialSupportDegs)
    instancesNum = size(emPosArcMin,1);
    timeBins = size(emPosArcMin,2);
    extraBins = 10;  % zero padding in time domain
    
    stimulusSequence     = zeros(instancesNum, timeBins+extraBins*2, size(stimulusImage,2), size(stimulusImage,3));
    powerSpectralDensity = zeros(instancesNum, timeBins+extraBins*2, size(stimulusImage,2), size(stimulusImage,3));
    
    pixelSizeDegs = spatialSupportDegs(2)-spatialSupportDegs(1);
    
    for instanceNo = 1:instancesNum
        % emPath for this instance
        theEMPosPathDegs = squeeze(emPosArcMin(instanceNo,:,:))/60;
        
        % Stimulus image for this instance
        theStimulusImage = squeeze(stimulusImage(instanceNo,:,:));
        
        % The zero contrast frames before and after the eye movement path
        theNullStimulus = theStimulusImage*0 + theStimulusImage(1,1);
        for tBin = 1:extraBins
            stimulusSequence(instanceNo, tBin,:,:)          = theNullStimulus;
            stimulusSequence(instanceNo, timeBins+tBin,:,:) = theNullStimulus;
        end
        
        % The frames during the eye movement path
        for tBin = 1:timeBins
            shiftAmountPixels = theEMPosPathDegs(tBin,:)/pixelSizeDegs;
            stimulusSequence(instanceNo, tBin+extraBins,:,:) = shiftImage(theStimulusImage, shiftAmountPixels);
        end
        
        % Compute 3D power spectral density
        powerSpectralDensity(instanceNo,:,:,:) = abs(fftshift(fftn(squeeze(stimulusSequence(instanceNo, :,:,:)))));
    end
    
    % mean over instances
    powerSpectralDensity = squeeze(mean(powerSpectralDensity,1));
    
    % Sum over the Y-frequency slice (vertical grating)
    powerSpectralDensityXT = squeeze(sum(powerSpectralDensity,2));
    powerSpectralDensityYT = squeeze(sum(powerSpectralDensity,3));
    powerSpectralDensityXY = squeeze(sum(powerSpectralDensity,1));
    
    % Make zero frequency entry nan for easier visualization
    yo = size(powerSpectralDensityXT,1)/2+1;
    xo = size(powerSpectralDensityXT,2)/2+1;
    powerSpectralDensityXT(yo,xo) = nan;
    
    yo = size(powerSpectralDensityYT,1)/2+1;
    xo = size(powerSpectralDensityYT,2)/2+1;
    powerSpectralDensityYT(yo,xo) = nan;
    
    yo = size(powerSpectralDensityXY,1)/2+1;
    xo = size(powerSpectralDensityXY,2)/2+1;
    powerSpectralDensityXY(yo,xo) = nan;
    
    N = size(powerSpectralDensityXY,1);
    N = N/2;
    sfMaxCPD = 1/(2*pixelSizeDegs);
    spatialFrequencySupport = ((-N):(N-1)) * sfMaxCPD * 1/N;
    
    N = size(powerSpectralDensityXT,1);
    N = N/2;
    tfMaxHz = 1/(2*(timeAxis(2)-timeAxis(1)));
    tfSupport = ((-N):(N-1))*tfMaxHz*1/N;
    
    figure(222); clf;
    subplot(1,3,1);
    imagesc(tfSupport, spatialFrequencySupport, powerSpectralDensityXT');
    axis 'square'
    xlabel('temporal frequency (Hz)');
    ylabel('spatial frequency, X (c/deg)');
    set(gca, 'YLim', [-30 30], 'XLim', [-100 100], 'FontSize', 14);
    
    
    subplot(1,3,2);
    imagesc(tfSupport, spatialFrequencySupport, powerSpectralDensityYT');
    axis 'square'
    xlabel('temporal frequency (Hz)');
    ylabel('spatial frequency, Y (c/deg)');
    set(gca, 'YLim', [-30 30], 'XLim', [-100 100], 'FontSize', 14);
    
    subplot(1,3,3);
    imagesc(spatialFrequencySupport, spatialFrequencySupport, powerSpectralDensityXY);
    axis 'square'
    xlabel('spatial frequency, X (c/deg)');
    ylabel('spatial frequency, Y (c/deg)');
    set(gca, 'XLim', [-30 30], 'YLim', [-30 30], 'FontSize', 14);
    colormap(gray)
    drawnow;
    pause
    
end

function shiftedImage = shiftImage(originalImage, shiftAmount)
    shiftedImage = imtranslate(originalImage,shiftAmount,'FillValues',originalImage(1,1));
    plotImages = false;
    if (plotImages)
        figure(999);
        subplot(1,2,1);
        imshow(originalImage); axis 'image';
        subplot(1,2,2);
        imshow(shiftedImage); axis 'image';
        drawnow
    end
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
    
    computePower = false;
    if (computePower)
        % Compute power
        gratingImagePower = computeImageTotalPower(gratingComponent);
    end
    
    % Generate noise component
    for instanceNo = 1:instancesNum
        noiseComponent = generateNoiseComponent(spatialSupportDegs, noiseParams);
 
        if (computePower)
            % Compute power
            noiseImagePower(instanceNo) = computeImageTotalPower(noiseComponent);
        end
        
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
    
    if (computePower)
        fprintf('Noise image power: %f\n', mean(noiseImagePower));
        fprintf('grating image power: %f\n', gratingImagePower);
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

function noiseImagePowerDB = computeImageTotalPower(img)
     noiseImagePowerDB = 10*log10(sum(sum((abs(fft2(img))).^2)));
end

function [gratingImage, gaussianEnvelope] = generateGratingComponent(spatialSupportDegs, gratingParams)

    N = length(spatialSupportDegs);
    [xx,yy] = meshgrid(spatialSupportDegs, spatialSupportDegs);
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
    
    visualizedInstancesNum = min([5 size(stimulusImage,1)]);
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

function [emPosArcMin, timeAxis] = generateFEM(nTrials, emDurationSeconds)
    sampleTimeSeconds = 1.0 / 1000;

    microSaccadeType = 'none'; % , 'heatmap/fixation based', 'none'
    fixEMobj = fixationalEM();
    fixEMobj.randomSeed = 1;
    fixEMobj.microSaccadeType = microSaccadeType;
    
    computeVelocity = ~true;
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, ...
            nTrials, computeVelocity, 'useParfor', true);
        
    timeAxis = fixEMobj.timeAxis;
    emPosArcMin = fixEMobj.emPosArcMin;
    
    % Center at (0,0)
    meanEMPos = squeeze(mean(emPosArcMin,2));
    for iTrial = 1:nTrials
        emPosArcMin(iTrial,:,1) = emPosArcMin(iTrial,:,1) - meanEMPos(iTrial,1);
        emPosArcMin(iTrial,:,2) = emPosArcMin(iTrial,:,2) - meanEMPos(iTrial,2);
    end
    
end

function plotFEM(emPosArcMin, timeAxis)
    trialNo = 1;
    xPosArcMin = squeeze(emPosArcMin(trialNo,:,1));
    yPosArcMin = squeeze(emPosArcMin(trialNo,:,2));
    
    maxPos = max(abs(emPosArcMin(:)));
    if (maxPos == 0)
        maxPos = 1;
    end
    posRange = maxPos*[-1 1];
    
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
