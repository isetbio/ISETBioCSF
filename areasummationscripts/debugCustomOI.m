function debugCustomOI

    subjectsNum = 256;
    calcPupilDiameterMM = 3.0;
    opticsModel = 'WvfHumanMeanOTFmagZeroOTFphase'; %'WvfHumanMeanOTFmagMeanOTFphase';
    theOI = generateMeanOTFcustomOI(opticsModel, calcPupilDiameterMM, subjectsNum);
    
    visualizeThePSFs(theOI);
    
end

function visualizeThePSFs(theOI)

    optics = oiGet(theOI,'optics');
    wavelengthsList = opticsGet(optics,'wave');
    fullOTF = opticsGet(optics, 'otf data');
    xSfGridCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
    ySfGridCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
    [xSfGridCyclesDegGrid,ySfGridCyclesDegGrid] = meshgrid(xSfGridCyclesDeg,ySfGridCyclesDeg);
    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 2538 1289]);
    
    for waveIndex = 1:numel(wavelengthsList)
        theWaveOTF = fullOTF(:,:,waveIndex);
        [xGridMinutes,yGridMinutes,theWavePSF] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(theWaveOTF));

        subplot(4,8,waveIndex);
        contourLevels = 0:0.01:1.0;
        contourf(xGridMinutes, yGridMinutes, theWavePSF/max(theWavePSF(:)), contourLevels, 'LineColor', 'none');
        hold on;

        xAx = xGridMinutes(1,:);
        yAx = yGridMinutes(:,1);
        centerPosition = floor(length(xAx)/2)+1;
        plot(xAx(centerPosition)*[1 1], [yAx(1) yAx(end)], 'r-');
        plot([xAx(1) xAx(end)], yAx(centerPosition)*[1 1], 'r-');
        axis 'square'
        title(sprintf('%d nm, vol: %1.9f',  wavelengthsList(waveIndex), sum(theWavePSF(:))));
        drawnow;
        colormap(jet(1024)); 
    end
    
end


function theOI = generateMeanOTFcustomOI(opticsModel, calcPupilDiameterMM, subjectsNum)

    theOI = oiCreate('wvf human', calcPupilDiameterMM);
    optics = oiGet(theOI,'optics');
    wavelengthsListToCompute = opticsGet(optics,'wave');
    xSfGridCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
    ySfGridCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
    [xSfGridCyclesDegGrid,ySfGridCyclesDegGrid] = meshgrid(xSfGridCyclesDeg,ySfGridCyclesDeg);
   
    alignmentWavelength = 550;
    [~,alignmentWIndex] = min(abs(wavelengthsListToCompute-alignmentWavelength));
    
    measPupilDiameterMM = 7.5;
    [Zcoeffs_SampleMean, Zcoeffs_S] = wvfLoadThibosVirtualEyes(measPupilDiameterMM);
    
    % Generate Z-coeffs for subjectsNum
    Zcoeffs = ieMvnrnd(Zcoeffs_SampleMean, Zcoeffs_S, subjectsNum)'; 
    
    umPerDegree = 300;
    tic
    for subjectIndex = 1:subjectsNum
        fprintf('Subject %d of %d\n', subjectIndex,subjectsNum);
        % Make WVF for this subject      
        wvfSubject = makeWVF(Zcoeffs(:,subjectIndex), wavelengthsListToCompute, ...
            measPupilDiameterMM, calcPupilDiameterMM, umPerDegree, sprintf('subject-%d', subjectIndex));
        
        % Compute translation vector for the PSF at the alignment wavelength
        translationVector = [];
        singleWavelengthOTF = wvfGet(wvfSubject, 'otf', wavelengthsListToCompute(alignmentWIndex)); 
        singleWavelengthPSF = wvfGet(wvfSubject, 'psf', wavelengthsListToCompute(alignmentWIndex));
        [~, translationVector] = otfWithZeroCenteredPSF(singleWavelengthOTF, singleWavelengthPSF, translationVector, xSfGridCyclesDegGrid,ySfGridCyclesDegGrid);
                
        for wIndex = 1:numel(wavelengthsListToCompute)
            singleWavelengthOTF = wvfGet(wvfSubject, 'otf', wavelengthsListToCompute(wIndex)); %opticsGet(optics, 'otf');
            singleWavelengthPSF = wvfGet(wvfSubject, 'psf', wavelengthsListToCompute(wIndex));
            
            % Center the OTF at this wavelength
            theCenteredOTF = otfWithZeroCenteredPSF(singleWavelengthOTF, singleWavelengthPSF, translationVector,  xSfGridCyclesDegGrid,ySfGridCyclesDegGrid);
            
            if (subjectIndex == 1) && (wIndex == 1)
                centeredOTFs = zeros(subjectsNum, size(theCenteredOTF,1), size(theCenteredOTF,2), numel(wavelengthsListToCompute));
            end
            centeredOTFs(subjectIndex, :,:, wIndex) = theCenteredOTF;
            
            visualize = false;
            if (visualize)
                figure(999);
                subplot(1,2,1);
                imagesc(fftshift(abs(singleWavelengthOTF)))
                axis 'square'
                subplot(1,2,2);
                imagesc(fftshift(abs(theCenteredOTF)))
                axis 'square'
                
                [xGridMinutes,yGridMinutes,psfFromOTF] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(singleWavelengthOTF));
                [~,~,centeredPSF] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(theCenteredOTF));

                figure(10000); clf;
                contourLevels = 0:0.1:1.0;

                subplot(2,2,1);
                contourf(xGridMinutes(1,:), xGridMinutes(1,:), singleWavelengthPSF/max(singleWavelengthPSF(:)), contourLevels, 'LineColor', 'none');
                axis 'square';
                title(sprintf('PSF direct from WVF (sum: %1.9f)', sum(singleWavelengthPSF(:))));

                subplot(2,2,2);
                contourf(xGridMinutes(1,:), xGridMinutes(1,:), psfFromOTF/max(psfFromOTF(:)), contourLevels, 'LineColor', 'none');
                axis 'square';
                title(sprintf('PSF indirect from OTF (sum: %1.9f)', sum(psfFromOTF(:))));

                subplot(2,2,3);
                contourf(xGridMinutes(1,:), xGridMinutes(1,:), centeredPSF/max(centeredPSF(:)), contourLevels, 'LineColor', 'none');
                axis 'square';
                title(sprintf('Centered PSF (sum: %1.9f)', sum(centeredPSF(:))));
        
                drawnow;
            end
            
        end % wIndex
    end % subjectIndex
    
    meanOTF = squeeze(mean(centeredOTFs,1));
    toc
    
    plotOTF = false;
    if (plotOTF)
        figure(1); clf;
    end
    
    for waveIndex = 1:numel(wavelengthsListToCompute)
        theWaveOTF = squeeze(meanOTF(:,:,waveIndex));
        
        otfMeanMag = abs(theWaveOTF);
        otfMeanPhase = angle(theWaveOTF);
        
        if strcmp(opticsModel, 'WvfHumanMeanOTFmagMeanOTFphase')
           phaseToUse = otfMeanPhase;
        elseif strcmp(opticsModel, 'WvfHumanMeanOTFmagZeroOTFphase')
           phaseToUse = 0*otfMeanPhase;
        end
        
        theWaveOTF = otfMeanMag .* exp(1i*phaseToUse);
        meanOTF(:,:,waveIndex) = theWaveOTF;
        
        if (plotOTF)
            [xGridMinutes,yGridMinutes,theWavePSF] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(theWaveOTF));

            subplot(4,8,waveIndex);
            contourLevels = 0:0.05:1.0;
            contourf(xGridMinutes, yGridMinutes, theWavePSF/max(theWavePSF(:)), contourLevels, 'LineColor', 'none');
            hold on;

            xAx = xGridMinutes(1,:);
            yAx = yGridMinutes(:,1);
            centerPosition = floor(length(xAx)/2)+1;
            plot(xAx(centerPosition)*[1 1], [yAx(1) yAx(end)], 'r-');
            plot([xAx(1) xAx(end)], yAx(centerPosition)*[1 1], 'r-');
            axis 'square'
            title(sprintf('%d nm, vol: %1.9f',  wavelengthsListToCompute(waveIndex), sum(theWavePSF(:))));
            colormap(gray(1024));
        end
    end
    

    % Update optics with new OTF data
    optics = opticsSet(optics,'otf data',meanOTF);

    % Stick optics back into oi
    theOI = oiSet(theOI,'optics',optics);
end

function theWVF = makeWVF(zcoeffs, wavelengthsToCompute, measPupilDiameterMM, calcPupilDiameterMM, umPerDegree, name)
    theWVF = wvfCreate(...
    			'umPerDegree', umPerDegree, ...
                'calc wavelengths',wavelengthsToCompute,...
                'zcoeffs',zcoeffs,...
                'name', name);
    
    % Set measured pupil size
    theWVF = wvfSet(theWVF,'measured pupil size', measPupilDiameterMM);
    
    % And pupil size for calculation
    theWVF = wvfSet(theWVF,'calc pupil size',calcPupilDiameterMM);
    
    % Now compute the PSF
    theWVF = wvfComputePSF(theWVF);
end

function [centeredOTF, translationVector] = otfWithZeroCenteredPSF(OTF, PSF, translationVector,  xSfGridCyclesDegGrid, ySfGridCyclesDegGrid)
    if (isempty(translationVector))
        % Compute center of mass
        centerOfMass = computeCenterOfMass(PSF);
        
        centerPosition = floor(size(PSF,1)/2) + 1 * [1 1];
        showTranslation = false;
        
        if (showTranslation)
            figure(200);clf
            imagesc(1:size(PSF,2), 1:size(PSF,1), PSF);
            hold on;
            plot(centerPosition(1)*[1 1], [1 size(PSF,1)], 'k-');
            plot([1 size(PSF,2)], centerPosition(2)*[1 1], 'k-');
            plot(centerOfMass(1), centerOfMass(2), 'ro', 'MarkerFaceColor', [1 0.5 0.5]);
            plot([centerPosition(1) centerOfMass(1)], [centerPosition(2)  centerOfMass(2)], 'r-');
            hold off;
            axis 'square'
            colormap(gray(1024));
            drawnow;
        end
        
        % Compute translation vector to bring center of mass at 0,0
        translationVector = (centerPosition-centerOfMass);
    end
    centeredOTF = shiftInFTplane(OTF, translationVector);

    % Make sure OTF(1,1) is equal to 1.0
    %[centeredOTF(1,1) max(abs(centeredOTF(:)))]
    %centeredOTF = centeredOTF/abs(centeredOTF(1,1));
    
    [xGridMinutes,yGridMinutes,centeredPSF] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(centeredOTF));
    
    if (abs(sum(centeredPSF(:))-1) > 1000*eps(1))
        fprintf('centeredPSF min = %1.9f, max = %1.9f, sum = %1.9f\n', min(centeredPSF(:)), max(centeredPSF(:)), sum(centeredPSF(:)));
        fprintf(2,'Sum does not equal to 1\n');
    end
   
end


function otf = shiftInFTplane(otf, translationVector)
    [N, M] = size(otf);
    x = [0:floor(N/2) floor(-N/2)+1:-1];
    y = [0:floor(M/2) floor(-M/2)+1:-1];
    x_shift = exp(-1i * 2 * pi * translationVector(2) * x' / N);
    y_shift = exp(-1i * 2 * pi * translationVector(1) * y / M);

    % Force conjugate symmetry. Otherwise this frequency component has no
    % corresponding negative frequency to cancel out its imaginary part.
    if mod(N, 2) == 0
        x_shift(N/2+1) = real(x_shift(N/2+1));
    end 
    if mod(M, 2) == 0
        y_shift(M/2+1) = real(y_shift(M/2+1));
    end
    otf = otf .* (x_shift * y_shift);
end

function centerOfMass = computeCenterOfMass(PSF)
    property = 'WeightedCentroid';
    c = regionprops(true(size(PSF)), PSF, property);
    centerOfMass = c.(property);
end