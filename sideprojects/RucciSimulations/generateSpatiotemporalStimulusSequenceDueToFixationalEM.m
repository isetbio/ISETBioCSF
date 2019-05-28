function sData = generateSpatiotemporalStimulusSequenceDueToFixationalEM(timeAxis, emPosArcMin, stimulus)
    instancesNum = size(emPosArcMin,1);
    timeBins = size(emPosArcMin,2);
    extraBins = 10;  % bins for zero padding in time domain
    totalTimeBins = timeBins+extraBins*2;
    
    % Preallocate memory
    sData.stimulusSequences = zeros(instancesNum, totalTimeBins, ...
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
        radius = floor(nRows/4);
        wSpatial = (welchWindow(2*radius))' * welchWindow(2*radius);

        
        nn = size(wSpatial,1);
        mm = size(wSpatial,2);
        oris = 0:45:315;
        timeSegmentsNum = 6;
        wTemporal = welchWindow(round(totalTimeBins/(timeSegmentsNum/2)));
        ll = numel(wTemporal);
        nEstimates = (length(oris)+1)*timeSegmentsNum
        xytWindow = zeros(nEstimates,totalTimeBins, size(stimulus.image,2), size(stimulus.image,3), 'single');
        
        % Overlapping spatiotemporal windowing for different spectral estimates
        
        k = 0;
        for timeSegment = 1:timeSegmentsNum
            wTemporal2 = zeros(1,totalTimeBins);
            wTemporal2((timeSegment-1)*floor(ll/2) + (1:ll)) = wTemporal;
            figure(444);
            subplot(timeSegmentsNum,1,timeSegment);
            plot(wTemporal2);
            drawnow;
            for orientation = 1:(numel(oris)+1)
                if (orientation < numel(oris)+1)
                    center(1) = nRows/2 + round(cosd(oris(orientation))*radius*0.8 - size(wSpatial,1)/2);
                    center(2) = mCols/2 + round(sind(oris(orientation))*radius*0.8 - size(wSpatial,2)/2);
                else
                    center = [nRows mCols]/2 - size(wSpatial)/2;
                end
                wSpatial2 = zeros(nRows, mCols);
                wSpatial2(center(1)+(1:size(wSpatial,1)), center(2)+(1:size(wSpatial,2))) = wSpatial;
                size(wSpatial2)
                figure(445);
                subplot(4,4,orientation);
                imagesc(wSpatial2);
                axis 'image';
                drawnow
                k = k + 1;
                for tBin = 1:totalTimeBins
                    xytWindow(k,tBin,:,:) = single(wSpatial2 * wTemporal2(tBin));
                end
            end
        end
        
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
        sData.stimulusSequences(instanceNo, :,:,:) = single(XYTstim);
        
        switch PSDmethod
            case 'FFT'
                powerSpectralDensity(instanceNo,:,:,:) = ...
                    1/(nRows*mCols*timeBins) * (abs(fftshift(fftn(XYTstim,2.^m)))).^2; 
                
            case 'windowedFFT'
                for k = 1:nEstimates
                    XYTstimTemp = XYTstim .* squeeze(double(xytWindow(k,:,:,:)));
                    if (k == 1)
                        spectralEstimate = 1/(nRows*mCols*timeBins) * (abs(fftshift(fftn(XYTstimTemp,2.^m)))).^2;
                    else
                        spectralEstimate = spectralEstimate + ...
                                        1/(nRows*mCols*timeBins) * (abs(fftshift(fftn(XYTstimTemp,2.^m)))).^2;
                    end
                end
                powerSpectralDensity(instanceNo,:,:,:) = spectralEstimate/nEstimates;
                
            otherwise
                error('Unknown PSD method: ''%s''. ', PSDmethod);
        end % switch
    end % for
    
    % mean spectral density over instances
    sData.meanSpatioTemporalSpectalDensity = squeeze(mean(powerSpectralDensity,1));
    
    N = size(sData.meanSpatioTemporalSpectalDensity,3);
    N = N/2;
    sfMaxCPD = 1/(2*pixelSizeDegs);
    sData.spatialFrequencySupport = ((-N):(N-1))/N * sfMaxCPD;
    
    N = size(sData.meanSpatioTemporalSpectalDensity,1);
    N = N/2;
    tfMaxHz = 1/(2*(timeAxis(2)-timeAxis(1)));
    sData.tfSupport = ((-N):(N-1))/N * tfMaxHz;
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