function stimulusSequence = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, stimulus, figNo)
    instancesNum = size(emPosArcMin,1);
    timeBins = size(emPosArcMin,2);
    extraBins = 50;  % bins for zero padding in time domain
    totalTimeBins = timeBins+extraBins*2;
    
    % Preallocate memory
    stimulusSequence = zeros(instancesNum, totalTimeBins, ...
        size(stimulus.image,2), size(stimulus.image,3), ...
        'single');

    pixelSizeDegs = stimulus.spatialSupportDegs(2)-stimulus.spatialSupportDegs(1);
    
    m = nextpow2([totalTimeBins size(stimulus.image,2), size(stimulus.image,3)]);
    powerSpectralDensity = zeros(instancesNum, 2^m(1), 2^m(2), 2^m(3));
                
    PSDmethod = 'FFT';
    PSDmethod = 'windowedFFT';
    
    if (strcmp(PSDmethod, 'windowedFFT'))
        nRows = size(stimulus.image,2);
        mCols = size(stimulus.image,3);
        wSpatial = (welchWindow(nRows))' * welchWindow(mCols);
        wTemporal = welchWindow(totalTimeBins);
        
        % Spatial windowing for each frame
        xytWindow = zeros(totalTimeBins, size(stimulus.image,2), size(stimulus.image,3));
        rectWindow = xytWindow + 1;
        for tBin = 1:totalTimeBins
            xytWindow(tBin,:,:) = wSpatial * wTemporal(tBin);
        end
        
        % Compute correction factor due to window
        winInt = sum(sum(sum(xytWindow.^2,1),2),3);
        rectInt = sum(sum(sum(rectWindow.^2,1),2),3);
        xytWindow = xytWindow  * rectInt/winInt;
    end
    
    for instanceNo = 1:instancesNum
        % emPath for this instance
        theEMPosPathDegs = squeeze(emPosArcMin(instanceNo,:,:))/60;
        
        % Stimulus image for this noise component
        theStimulusImage = squeeze(stimulus.image(instanceNo,:,:));
        
        % The zero contrast frames before and after the eye movement path
        theNullStimulus = theStimulusImage*0 + theStimulusImage(1,1);
        
        % Pleallocate memory
        XYTstim = zeros(totalTimeBins, size(theNullStimulus,1), size(theNullStimulus,2));
        % The null frames  before and after the eye movement path
        for tBin = 1:extraBins
            XYTstim(tBin,:,:)          = theNullStimulus;
            XYTstim(timeBins+tBin,:,:) = theNullStimulus;
        end
        
        % The frames during the eye movement path
        for tBin = 1:timeBins
            shiftAmountPixels = theEMPosPathDegs(tBin,:)/pixelSizeDegs;
            XYTstim(extraBins+tBin,:,:) = shiftImage(theStimulusImage, shiftAmountPixels);
        end
        
        % spatiotemporal movie of this noise instance
        stimulusSequence(instanceNo, :,:,:) = single(XYTstim);
        
        switch PSDmethod
            case 'FFT'
                powerSpectralDensity(instanceNo,:,:,:) = ...
                    1/(nRows*mCols) * (abs(fftshift(fftn(XYTstim,2.^m)))).^2; 
                
            case 'windowedFFT'
                XYTstim = XYTstim .* xytWindow;
                powerSpectralDensity(instanceNo,:,:,:) = ...
                    1/(nRows*mCols) * (abs(fftshift(fftn(XYTstim,2.^m)))).^2; 
                
            otherwise
                error('Unknown PSD method: ''%s''. ', PSDmethod);
        end % switch
    end % for
    
    % mean spectral density over instances
    powerSpectralDensity = squeeze(mean(powerSpectralDensity,1));
    
    [~,idx] = max(powerSpectralDensity(:));
    [maxTFindex, maxSFYindex, maxSFXindex] = ind2sub(size(powerSpectralDensity),idx);
    
    % Mean over the temporal frequency axis
    powerSpectralDensityXY = squeeze(mean(powerSpectralDensity,1));
    powerSpectralDensityXYpeakTF = squeeze(powerSpectralDensity(maxTFindex,:,:));
    
    % Mean over the Y-spatial frequency axis
    powerSpectralDensityXT = squeeze(mean(powerSpectralDensity,2));
    powerSpectralDensityXTpeakSFY = squeeze(powerSpectralDensity(:,maxSFYindex,:));
    
    % Mean over the X-spatial frequency axis
    powerSpectralDensityYT = squeeze(mean(powerSpectralDensity,3));
    powerSpectralDensityYTpeakSFX = squeeze(powerSpectralDensity(:,:,maxSFXindex));
    
    
    % Make zero frequency entry nan for easier visualization
%     yo = size(powerSpectralDensityXT,1)/2+1;
%     xo = size(powerSpectralDensityXT,2)/2+1;
%     powerSpectralDensityXT(yo,xo) = nan;
%     
%     yo = size(powerSpectralDensityYT,1)/2+1;
%     xo = size(powerSpectralDensityYT,2)/2+1;
%     powerSpectralDensityYT(yo,xo) = nan;
%     
%     yo = size(powerSpectralDensityXY,1)/2+1;
%     xo = size(powerSpectralDensityXY,2)/2+1;
%     powerSpectralDensityXY(yo,xo) = nan;
    
    N = size(powerSpectralDensityXY,1);
    N = N/2;
    sfMaxCPD = 1/(2*pixelSizeDegs);
    spatialFrequencySupport = ((-N):(N-1))/N * sfMaxCPD;
    
    N = size(powerSpectralDensityXT,1);
    N = N/2;
    tfMaxHz = 1/(2*(timeAxis(2)-timeAxis(1)));
    tfSupport = ((-N):(N-1))/N * tfMaxHz;
    
    fprintf('Max PSD at %2.1f Hz, %2.1f c/deg (x), %2.1f c/deg (y)\n', ...
        tfSupport(maxTFindex), spatialFrequencySupport(maxSFXindex), spatialFrequencySupport(maxSFYindex));
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1100 320]);
    
    sfLims = [0 15];
    tfLims = [-30 30];
    logScaling = true;
    
    if (logScaling)
        visualizedRangeDeciBels = 40;
        powerSpectralDensityXT = 10*log10(powerSpectralDensityXT);
        powerSpectralDensityYT = 10*log10(powerSpectralDensityYT);
        powerSpectralDensityXY = 10*log10(powerSpectralDensityXY);
        powerSpectralDensityYTpeakSFX = 10*log10(powerSpectralDensityYTpeakSFX);
        powerSpectralDensityXTpeakSFY = 10*log10(powerSpectralDensityXTpeakSFY);
        powerSpectralDensityXYpeakTF = 10*log10(powerSpectralDensityXYpeakTF);
        cLims = [-visualizedRangeDeciBels 0] + max([max(powerSpectralDensityXT) max(powerSpectralDensityYT) max(powerSpectralDensityXY)]);
        cLims2 = [-visualizedRangeDeciBels 0] + max(10*log10(powerSpectralDensity(:)));
    else
        cLims = [0 max(powerSpectralDensity(:))];
        cLims2 = cLims;
    end
    
    subplot(2,3,1);
    imagesc(tfSupport, spatialFrequencySupport, powerSpectralDensityXT'); hold on;
    plot([0 0], [spatialFrequencySupport(1) spatialFrequencySupport(end)], 'g-');
    plot([tfSupport(1) tfSupport(end)], [0 0], 'r-');
    axis 'square'
    colorbar('northoutside')
    title('mean PSD over Y-spatial frequency');
    xlabel('temporal frequency (Hz)');
    ylabel('spatial frequency, X (c/deg)');
    set(gca, 'XLim', tfLims, 'YLim', sfLims ,  'FontSize', 14);
    set(gca, 'CLim', cLims);
    
    subplot(2,3,2);
    imagesc(tfSupport, spatialFrequencySupport, powerSpectralDensityYT'); hold on;
    plot([0 0], [spatialFrequencySupport(1) spatialFrequencySupport(end)], 'g-');
    plot([tfSupport(1) tfSupport(end)], [0 0], 'r-');
    axis 'square'
    colorbar('northoutside')
    title('mean PSD over X-spatial frequency');
    xlabel('temporal frequency (Hz)');
    ylabel('spatial frequency, Y (c/deg)');
    set(gca, 'XLim', tfLims, 'YLim', sfLims , 'FontSize', 14);
    set(gca, 'CLim', cLims);
    
    subplot(2,3,3);
    imagesc(spatialFrequencySupport, spatialFrequencySupport, powerSpectralDensityXY);
    hold on
    plot([0 0], [spatialFrequencySupport(1) spatialFrequencySupport(end)], 'g-');
    plot([spatialFrequencySupport(1) spatialFrequencySupport(end)], [0 0], 'g-');
    axis 'square'
    colorbar('northoutside')
    title('mean PSD over temporal frequency');
    xlabel('spatial frequency, X (c/deg)');
    ylabel('spatial frequency, Y (c/deg)');
    set(gca, 'XLim', sfLims(2)*[-1 1], 'YLim', sfLims, 'FontSize', 14);
    set(gca, 'CLim', cLims);
    
    
    subplot(2,3,4);
    imagesc(tfSupport, spatialFrequencySupport, powerSpectralDensityXTpeakSFY'); hold on;
    plot([0 0], [spatialFrequencySupport(1) spatialFrequencySupport(end)], 'g-');
    plot([tfSupport(1) tfSupport(end)], [0 0], 'r-');
    axis 'square'
    colorbar('northoutside')
    title(sprintf('at sfY = %2.1f c/deg',spatialFrequencySupport(maxSFYindex)));
    xlabel('temporal frequency (Hz)');
    ylabel('spatial frequency, X (c/deg)');
    set(gca, 'XLim', tfLims, 'YLim', sfLims ,  'FontSize', 14);
    set(gca, 'CLim', cLims2);
    
    
    subplot(2,3,5);
    imagesc(tfSupport, spatialFrequencySupport, powerSpectralDensityYTpeakSFX'); hold on;
    plot([0 0], [spatialFrequencySupport(1) spatialFrequencySupport(end)], 'g-');
    plot([tfSupport(1) tfSupport(end)], [0 0], 'r-');
    axis 'square'
    colorbar('northoutside')
    title(sprintf('at sfX = %2.1f c/deg',spatialFrequencySupport(maxSFXindex)));
    xlabel('temporal frequency (Hz)');
    ylabel('spatial frequency, Y (c/deg)');
    set(gca, 'XLim', tfLims, 'YLim', sfLims , 'FontSize', 14);
    set(gca, 'CLim', cLims2);
    
    subplot(2,3,6);
    imagesc(spatialFrequencySupport, spatialFrequencySupport, powerSpectralDensityXYpeakTF);
    hold on
    plot([0 0], [spatialFrequencySupport(1) spatialFrequencySupport(end)], 'g-');
    plot([spatialFrequencySupport(1) spatialFrequencySupport(end)], [0 0], 'g-');
    axis 'square'
    colorbar('northoutside')
    title(sprintf('at tf = %2.1f Hz',tfSupport(maxTFindex)));
    xlabel('spatial frequency, X (c/deg)');
    ylabel('spatial frequency, Y (c/deg)');
    set(gca, 'XLim', sfLims(2)*[-1 1], 'YLim', sfLims, 'FontSize', 14);
    set(gca, 'CLim', cLims2);
    
    colormap(hot)
    drawnow;
    
end

function w = welchWindow(N)
     w  = (1 - (((0:N-1)-((N-1)/2))/((N+1)/2)).^2);
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


