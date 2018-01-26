function [theCustomOI, Zcoeffs] = oiWithCustomOptics(opticsModel, wavefrontSpatialSamples, calcPupilDiameterMM, umPerDegree)

    % Save rng state. 
    currentState = rng;
    
    % Set rng seed to 1 to ensure we always get the same Zcoeffs
    rng(1);

    switch opticsModel
        case  'AOoptics75mmPupil'
            subjectsNum = 1;
            subjectIndices = 1;
            calcPupilDiameterMM = 7.5;
            
        case {'WvfHumanMeanOTFmagMeanOTFphase', 'WvfHumanMeanOTFmagZeroOTFphase'}
            subjectsNum = 256;
            subjectIndices = 1:subjectsNum;

        case 'WvfHumanSubject1'
             subjectsNum = 100;
             subjectIndices = 1;

        case 'WvfHumanSubject2'
             subjectsNum = 100;
             subjectIndices = 17;
             
        case 'WvfHumanSubject3'
             subjectsNum = 100;
             subjectIndices = 55;
             
        case 'WvfHumanSubject4'
             subjectsNum = 100;
             subjectIndices = 69;
             
        case 'WvfHumanSubject5'
             subjectsNum = 100;
             subjectIndices = 78;     
    end
    
    showTranslation = false;
    centerPSF = true;
    [theCustomOI, Zcoeffs] = generateMeanAcrossSubjectsOTFcustomOI(opticsModel, wavefrontSpatialSamples, calcPupilDiameterMM, subjectsNum, subjectIndices, umPerDegree, centerPSF, showTranslation);
    
    % Restore current rng state
    rng(currentState);
end


function [theCustomOI, Zcoeffs] = generateMeanAcrossSubjectsOTFcustomOI(opticsModel, wavefrontSpatialSamples, calcPupilDiameterMM, subjectsNum, subjectIndices, umPerDegree, centerPSF, showTranslation)

    theOI = oiCreate('wvf human', calcPupilDiameterMM);
    optics = oiGet(theOI,'optics');
    wavelengthsListToCompute = opticsGet(optics,'wave');
   
    % Alignment wavelength
    alignmentWavelength = 550;
    [~,alignmentWIndex] = min(abs(wavelengthsListToCompute-alignmentWavelength));
    
    % Load the Thibos data set.  If we load 7.5 mm, we can compute for anything possible.
    measPupilDiameterMM = 7.5;
    [Zcoeffs_SampleMean, Zcoeffs_S] = wvfLoadThibosVirtualEyes(measPupilDiameterMM);
    
    Zcoeffs = ieMvnrnd(Zcoeffs_SampleMean, Zcoeffs_S, subjectsNum)'; 
    Zcoeffs = Zcoeffs(:,subjectIndices);
        
    if contains(opticsModel, 'AOoptics')
        Zcoeffs = Zcoeffs * 0;
    end

    tic
    for subjectIndex = 1:numel(subjectIndices)
        fprintf('Generating wavefront object for subject %d of %d\n', subjectIndex,numel(subjectIndices));
        
        % Make a WVF object for this subject      
        wvfSubject = makeWVF(wavefrontSpatialSamples, Zcoeffs(:,subjectIndex), wavelengthsListToCompute, ...
            measPupilDiameterMM, calcPupilDiameterMM, umPerDegree, sprintf('subject-%d', subjectIndex));
        
        % Compute an optics structure from the wvf object
        optics = oiGet(wvf2oi(wvfSubject),'optics');
        
        if contains(opticsModel, 'AOoptics')
            theCustomOI = oiSet(theOI,'optics', optics);
            continue;
        end
        
        % Proceed with non AOoptics
        xSfGridCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
        ySfGridCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
        [xSfGridCyclesDegGrid,ySfGridCyclesDegGrid] = meshgrid(xSfGridCyclesDeg,ySfGridCyclesDeg);

        allWavesOTF = opticsGet(optics, 'otf');
        allWavePSFs = wvfGet(wvfSubject, 'psf');

        if (centerPSF)
            % Compute translation vector for the PSF at the alignment wavelength
            translationVector = [];
            %singleWavelengthOTF = allWavesOTF(:,:,alignmentWIndex); 
            %singleWavelengthPSF = allWavePSFs{alignmentWindex};
            alignmentWavelengthOTF = wvfGet(wvfSubject, 'otf', wavelengthsListToCompute(alignmentWIndex)); %opticsGet(optics, 'otf');
            alignmentWavelengthPSF = wvfGet(wvfSubject, 'psf', wavelengthsListToCompute(alignmentWIndex));

            [~, translationVector, ~, ~, ~] = otfWithZeroCenteredPSF(alignmentWavelengthOTF, alignmentWavelengthPSF, translationVector, xSfGridCyclesDegGrid,ySfGridCyclesDegGrid, showTranslation);
        else
            fprintf(2,'Not centering PSF\n');
            translationVector = [0 0];
        end
        
        for wIndex = 1:numel(wavelengthsListToCompute)
            
            %singleWavelengthOTF = allWavesOTF(:,:,wIndex); 
            %singleWavelengthPSF = allWavePSFs{wIndex};
            
            theWaveOTFbeforeCentering = wvfGet(wvfSubject, 'otf', wavelengthsListToCompute(wIndex)); %opticsGet(optics, 'otf');
            theWavePSFbeforeCentering = wvfGet(wvfSubject, 'psf', wavelengthsListToCompute(wIndex));
            
            if (centerPSF)
                % Center the OTF at this wavelength
                [theWaveOTF, ~, theWavePSF, xGridMinutes,yGridMinutes] = otfWithZeroCenteredPSF(theWaveOTFbeforeCentering, theWavePSFbeforeCentering, translationVector,  xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showTranslation);
            else
                theWaveOTF = theWaveOTFbeforeCentering;
                [xGridMinutes,yGridMinutes,theWavePSF] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(theWaveOTF));
                theWavePSF = fftshift(ifft2(theWaveOTF)); % 
            end
            
            % Make sure OTF(1,1) is equal to 1.0
            [theWaveOTF(1,1) max(abs(theWaveOTF(:)))]
            theWaveOTF = theWaveOTF/abs(theWaveOTF(1,1));
    
            if (subjectIndex == 1) && (wIndex == 1)
                centeredOTFs = zeros(numel(subjectIndices), size(theWaveOTF,1), size(theWaveOTF,2), numel(wavelengthsListToCompute));
            end
            
            centeredOTFs(subjectIndex, :,:, wIndex) = theWaveOTF;
            
            visualizeDebug = true;
            if (visualizeDebug)
                figure(10000);
                subplot(2,2,1);
                imagesc(squeeze(xSfGridCyclesDegGrid(1,:)), squeeze(ySfGridCyclesDegGrid(:,1)), fftshift(abs(theWaveOTFbeforeCentering)));
                axis 'square'
                title(sprintf('OTF (%2.0f nm) mag before centering', wavelengthsListToCompute(wIndex)));
                
                subplot(2,2,2);
                imagesc(squeeze(xSfGridCyclesDegGrid(1,:)), squeeze(ySfGridCyclesDegGrid(:,1)), fftshift(abs(theWaveOTF)))
                axis 'square'
                title(sprintf('OTF (%2.0f nm) mag after centering', wavelengthsListToCompute(wIndex)));
                
                subplot(2,2,3);
                imagesc(xGridMinutes(1,:), yGridMinutes(:,1), theWavePSFbeforeCentering/max(theWavePSFbeforeCentering(:)), [0 1]);
                axis 'square';
                title(sprintf('PSF (%2.0f nm) before centering (sum: %1.9f)', wavelengthsListToCompute(wIndex), sum(theWavePSFbeforeCentering(:))));

                subplot(2,2,4);
                imagesc(xGridMinutes(1,:), yGridMinutes(:,1), theWavePSF/max(theWavePSF(:)), [0 1]);
                axis 'square';
                title(sprintf('PSF (%2.0f nm) after centering (sum: %1.9f)', wavelengthsListToCompute(wIndex), sum(theWavePSF(:)))); 
                drawnow;
                pause
            end
            
        end % wIndex
    end % subjectIndex
    
    if (contains(opticsModel, 'AOoptics'))
        fprintf('\nCompute of ''%s'' optics took %2.3f minutes.\n', opticsModel, toc/60);
        return;
    end
    
    % Mean across subjects
    meanOTF = squeeze(mean(centeredOTFs,1));

    
    % Reconstruct from mag and phase
    for waveIndex = 1:numel(wavelengthsListToCompute)
        theWaveOTF = squeeze(meanOTF(:,:,waveIndex));

        otfMeanMag = abs(theWaveOTF);
        otfMeanPhase = angle(theWaveOTF);

        phaseToUse = otfMeanPhase;
        if strcmp(opticsModel, 'WvfHumanMeanOTFmagZeroOTFphase')
           phaseToUse = 0*otfMeanPhase;
        end

        theWaveOTF = otfMeanMag .* exp(1i*phaseToUse);
        meanOTF(:,:,waveIndex) = theWaveOTF;
    end
    
    
    % Make new optics with new OTF data
    customOptics = opticsSet(optics,'otf data',meanOTF);

    % Create oi with custom optics 
    theCustomOI = oiSet(theOI,'optics', customOptics);
    
    fprintf('\nCompute of ''%s'' optics took %2.3f minutes.\n', opticsModel, toc/60);
    
end

function theWVF = makeWVF(wavefrontSpatialSamples, zcoeffs, wavelengthsToCompute, measPupilDiameterMM, calcPupilDiameterMM, umPerDegree, name)
    theWVF = wvfCreate(...
    			'umPerDegree', umPerDegree, ...
                'calc wavelengths',wavelengthsToCompute,...
                'measuredpupil', measPupilDiameterMM, ...
                'spatialsamples', wavefrontSpatialSamples, ...
                'zcoeffs', zcoeffs,...
                'name', name);
    
    theWVF = wvfCreate('calc wavelengths', [400:10:700]);
    % Set measured pupil size
    %theWVF = wvfSet(theWVF,'measured pupil size', measPupilDiameterMM);
    
    % And pupil size for calculation
    %theWVF = wvfSet(theWVF,'calc pupil size',calcPupilDiameterMM);
    
    % Set a high number of samples
    %theWVF = wvfSet(theWVF, 'number spatial samples', wavefrontSpatialSamples);
    
    % Now compute the PSF
    theWVF = wvfComputePSF(theWVF);
end

function [centeredOTF,  translationVector, centeredPSF, xGridMinutes,yGridMinutes] = otfWithZeroCenteredPSF(OTF, PSF, translationVector,  xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showTranslation)
    if (isempty(translationVector))
        % Compute center of mass
        centerOfMass = computeCenterOfMass(PSF);
        centerOfMassNative = computeCenterOfMassNative(PSF);
        centerPosition = floor(size(PSF,1)/2) + 1 * [1 1];
        fprintf('Center of mass: method1 (%f %f), method2 (%f,%f)\n', centerOfMass(1), centerOfMass(2), centerOfMassNative(1), centerOfMassNative(2));

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

    
    [xGridMinutes,yGridMinutes,centeredPSF] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(centeredOTF));
    
    if (abs(sum(centeredPSF(:))-1) > 1000*eps(1))
        fprintf('centeredPSF min = %1.9f, max = %1.9f, sum = %1.9f\n', min(centeredPSF(:)), max(centeredPSF(:)), sum(centeredPSF(:)));
        error('PSF volume does not equal 1\n');
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

function centerOfMass = computeCenterOfMassNative(PSF)
    [rc,cc] = ndgrid(1:size(PSF,1),1:size(PSF,2));
    Mt = sum(PSF(:));
    centerOfMassY = sum(PSF(:) .* rc(:)) / Mt;
    centerOfMassX = sum(PSF(:) .* cc(:)) / Mt;
    centerOfMass = [centerOfMassX centerOfMassY];
end
