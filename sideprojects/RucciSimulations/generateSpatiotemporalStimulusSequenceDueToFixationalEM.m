function stimulusSequence = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, stimulus, figNo)
    instancesNum = size(emPosArcMin,1);
    timeBins = size(emPosArcMin,2);
    extraBins = 0;  % zero padding in time domain
    
    stimulusSequence = zeros(instancesNum, timeBins+extraBins*2, size(stimulus.image,2), size(stimulus.image,3));

    pixelSizeDegs = stimulus.spatialSupportDegs(2)-stimulus.spatialSupportDegs(1);
    
    PSDmethod = 'FFT';
    PSDmethod = 'windowedFFT';
    
    for instanceNo = 1:instancesNum
        % emPath for this instance
        theEMPosPathDegs = squeeze(emPosArcMin(instanceNo,:,:))/60;
        
        % Stimulus image for this instance
        theStimulusImage = squeeze(stimulus.image(instanceNo,:,:));
        
        % The zero contrast frames before and after the eye movement path
        theNullStimulus = theStimulusImage*0 + theStimulusImage(1,1);
        for tBin = 1:extraBins
            stimulusSequence(instanceNo, tBin,:,:)          = theNullStimulus;
            stimulusSequence(instanceNo, timeBins+tBin,:,:) = theNullStimulus;
        end
        
        % The frames during the eye movement path
        for tBin = 1:timeBins
            shiftAmountPixels = theEMPosPathDegs(tBin,:)/pixelSizeDegs;
            stimulusSequence(instanceNo, extraBins+tBin,:,:) = shiftImage(theStimulusImage, shiftAmountPixels);
        end
        
        % spatiotemporal movie of this noise instance
        XYTstim = squeeze(stimulusSequence(instanceNo, :,:,:));
        
        switch PSDmethod
            case 'FFT'
                % Compute 3D power spectral density with no windowing (FFT)
                m = nextpow2(size(XYTstim));
                if (instanceNo == 1)
                    powerSpectralDensity = zeros(instancesNum, 2^m(1), 2^m(2), 2^m(3));
                end
                powerSpectralDensity(instanceNo,:,:,:) = (abs(fftshift(fftn(XYTstim,2.^m)))).^2;
            
            case 'windowedFFT'
                k = size(XYTstim,1);
                n = size(XYTstim,2);
                m = size(XYTstim,3);
                wSpatial = (1 - (((0:n-1)-((n-1)/2))/((n+1)/2)).^2)' * ...
                      (1 - (((0:m-1)-((m-1)/2))/((m+1)/2)).^2); % Welch window
                
                wTemporal = (1 - (((0:k-1)-((k-1)/2))/((k+1)/2)).^2);
                xytWindow = XYTstim * 0;
                
                % Spatial windowing for each frame
                for tBin = 1:size(XYTstim,1)
                    xytWindow(tBin,:,:) = wSpatial * wTemporal(tBin);
                end
                
                
                winInt = sum(sum(sum(xytWindow.^2,1),2),3);
                rectInt = prod(size(XYTstim));
                
                % Window
                XYTstim = XYTstim .* xytWindow;
                
                m = nextpow2(size(XYTstim));
                if (instanceNo == 1)
                    powerSpectralDensity = zeros(instancesNum, 2^m(1), 2^m(2), 2^m(3));
                end
                
                powerSpectralDensity(instanceNo,:,:,:) = 1/winInt * (abs(fftshift(fftn(XYTstim,2.^m)))).^2;
                
                
            otherwise
                error('Unknown PSD method: ''%s''. ', PSDmethod);
        end
    end
    
    % mean over instances
    powerSpectralDensity = squeeze(mean(powerSpectralDensity,1));
    
    % Mean over the Y-spatial frequency axis
    powerSpectralDensityXT = squeeze(mean(powerSpectralDensity,2));
    
    % Mean over the X-spatial frequency axis
    powerSpectralDensityYT = squeeze(mean(powerSpectralDensity,3));
    
    % Mean over the temporal frequency axis
    powerSpectralDensityXY = squeeze(mean(powerSpectralDensity,1));
    
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
    spatialFrequencySupport = ((-N):(N-1)) * sfMaxCPD * 1/N;
    
    N = size(powerSpectralDensityXT,1);
    N = N/2;
    tfMaxHz = 1/(2*(timeAxis(2)-timeAxis(1)));
    tfSupport = ((-N):(N-1))*tfMaxHz*1/N;
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1100 320]);
    
    sfLims = [0 15];
    tfLims = [-30 30];
    logScaling = true;
    
    if (logScaling)
        powerSpectralDensityXT = 10*log10(powerSpectralDensityXT);
        powerSpectralDensityYT = 10*log10(powerSpectralDensityYT);
        powerSpectralDensityXY = 10*log10(powerSpectralDensityXY);
        cLims = [-30 0] + max([max(powerSpectralDensityXT) max(powerSpectralDensityYT) max(powerSpectralDensityXY)]);
    else
        cLims = [0 max(powerSpectralDensity(:))];
    end
    
    subplot(1,3,1);
    
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
    
    subplot(1,3,2);
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
    
    subplot(1,3,3);
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
    
    colormap(hot)
    drawnow;
    
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


