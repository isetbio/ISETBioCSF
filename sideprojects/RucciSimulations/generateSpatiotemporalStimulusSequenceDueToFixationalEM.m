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
    
    fftSize = 2.^(nextpow2([totalTimeBins size(stimulus.image,2), size(stimulus.image,3)]));
    powerSpectralDensity = zeros(instancesNum, fftSize(1), fftSize(2), fftSize(3));
                
    PSDmethod = 'FFT';
    PSDmethod = 'windowedFFT';   
    
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
                    1/prod(fftSize) * (abs(fftshift(fftn(XYTstim,fftSize)))).^2; 
                
            case 'windowedFFT'
                nRows = size(stimulus.image,2);
                mCols = size(stimulus.image,3);
        
                radius = floor(nRows/4);
                spatialSamples = 2*radius;

                timeSegmentsNum = 6;
                temporalSamples = round(totalTimeBins/(timeSegmentsNum/2));

                oris = 0:45:315;
                
                powerSpectralDensity(instanceNo,:,:,:) = ...
                    welchSpectrum(XYTstim, timeSegmentsNum, temporalSamples, spatialSamples, oris, totalTimeBins, nRows, mCols, fftSize);
                
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

function powerSpectralDensity = welchSpectrum(XYTstim, timeSegmentsNum, temporalSamples, spatialSamples, oris, totalTimeBins, nRows, mCols, fftSize)

    radius = spatialSamples/2;
    nEstimates = 0;
    
    for timeSegment = 1:timeSegmentsNum
        temporalOffset = (timeSegment-1)*floor(temporalSamples/2);
            
        for orientation = 1:(numel(oris)+1)
            fprintf('Computing spectral estimate %d/%d,  %d/%d\n', timeSegment, timeSegmentsNum, orientation, (numel(oris)+1));
            
            if (orientation < numel(oris)+1)
                spatialOffset(1) = (nRows-spatialSamples)/2 + round(cosd(oris(orientation))*radius*0.8);
                spatialOffset(2) = (mCols-spatialSamples)/2 + round(sind(oris(orientation))*radius*0.8);
            else
                spatialOffset = [nRows-spatialSamples mCols-spatialSamples]/2;
            end
            
            [xytWelchWindow, wTemporal2, wSpatial2] = ...
                spatioTemporalWelchWindow(temporalSamples, spatialSamples, temporalOffset, spatialOffset, totalTimeBins, nRows, mCols);
            
            if (orientation == 1)
                figure(444);
                subplot(timeSegmentsNum,1,timeSegment);
                plot(wTemporal2);
                drawnow;
            end
            
            if (timeSegment == 1)
                figure(445);
                subplot(4,4,orientation);
                imagesc(wSpatial2);
                axis 'image';
                drawnow
            end
            
            % FFT of windowed spatiotemporal response
            thisSpectralEstimate = 1/prod(fftSize) * (abs(fftshift(fftn((XYTstim .* xytWelchWindow),fftSize)))).^2;
            
            nEstimates = nEstimates + 1;
            if (nEstimates == 1)
                spectralEstimate = thisSpectralEstimate;
            else
                spectralEstimate = spectralEstimate + thisSpectralEstimate;
            end
            
        end % orientation
    end % timeSegment
        
    powerSpectralDensity = spectralEstimate / nEstimates;
end

function [xytWindow, wTemporal2, wSpatial2] = ...
    spatioTemporalWelchWindow(temporalSamples, spatialSamples, temporalOffset, spatialOffset, tBins, nRows, mCols)
%
    wTemporal = welchWindow(temporalSamples);
    wTemporal2 = zeros(1,tBins);
    wTemporal2(temporalOffset + (1:length(wTemporal))) = wTemporal;
%
    wSpatial = (welchWindow(spatialSamples))' * welchWindow(spatialSamples);
    wSpatial2 = zeros(nRows, mCols);
    wSpatial2(spatialOffset(1)+(1:size(wSpatial,1)), spatialOffset(2)+(1:size(wSpatial,2))) = wSpatial;
    xytWindow = zeros(tBins, nRows, mCols);
%
    for tBin = 1:tBins
        xytWindow(tBin,:,:) = wSpatial2 * wTemporal2(tBin);
    end
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