function [theOI, theCustomOI] = oiWithCustomOptics(pupilDiameterMM, opticsModel)
    % 
    theOI = oiCreate('wvf human', pupilDiameterMM);
    theCustomOI = theOI;

    if (~ismember(opticsModel, {'WvfHumanMeanOTFmagMeanOTFphase', 'WvfHumanMeanOTFmagZeroOTFphase'}))
        error('Unknown optics model\n');
    end

    % Load the 4.5 mm pupil diameter Thibos data set
    measPupilDiameterMM = 4.5;
    [Zcoeffs_SampleMean, Zcoeffs_S] = wvfLoadThibosVirtualEyes(measPupilDiameterMM);

    subjectsNum  = 100;
    % Generate Z-coeffs for subjectsNum
    Zcoeffs = ieMvnrnd(Zcoeffs_SampleMean, Zcoeffs_S, subjectsNum)'; 

    % Pull out the optics structure
    optics = oiGet(theCustomOI,'optics');
    wavelengthsListToCompute = opticsGet(optics,'wave');

    customMTFdata = computeCustomMTF(Zcoeffs, ...
           pupilDiameterMM, wavelengthsListToCompute, subjectsNum, opticsModel);
    
    % Update optics with new OTF data
    optics = opticsSet(optics,'otf data',customMTFdata.otf);

    % Stick optics into oi
    theCustomOI = oiSet(theCustomOI,'optics',optics);
end

function customMTFdata = computeCustomMTF(Zcoeffs, pupilDiameterMM, wavelengths, subjectsNum, opticsModel)

    for subjectIndex = 1:subjectsNum
        % Make WVF for this subject
        
        wvfSubject = makeWVF(Zcoeffs(:,subjectIndex), wavelengths, pupilDiameterMM, sprintf('subject-%d', subjectIndex));
        optics = oiGet(wvf2oi(wvfSubject),'optics'); 
        otf = opticsGet(optics, 'otf');
        
        % Get PSF at selected wavelength
        referencePSF = wvfGet(wvfSubject, 'psf', 550);
        
        for waveIndex = 1:numel(wavelengths)
            % Get OTF at  wavelength
            theWaveOTF = squeeze(otf(:,:,waveIndex));
            
            if (waveIndex == 1)
                [theWaveOTF, translationVector] = otfWithZeroCenteredPSF(theWaveOTF, referencePSF, []);
            else
                theWaveOTF = otfWithZeroCenteredPSF(theWaveOTF, referencePSF, translationVector);
            end
            otf(:,:,waveIndex) = theWaveOTF;
            
            if (subjectIndex == 1) && (waveIndex == 1)
                otfSubjects = zeros(subjectsNum, size(otf,1), size(otf,2), size(otf,3));
                customMTFdata.xSfGridCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
                customMTFdata.ySfGridCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
                customMTFdata.otf = 0*otf;
            end
        end % waveIndex
        otfSubjects(subjectIndex,:,:,:) = otf;
    end % subjectIndex
    
    % Compute average OTF across subjects
    otfMean = squeeze(mean(otfSubjects,1));
    
    % Compute radially-symmetric OTF with desired phase at all wavelengths
    for waveIndex = 1:numel(wavelengths)
        theWaveOTF = squeeze(otfMean(:,:,waveIndex));
        otfMeanMag = fftshift(abs(theWaveOTF));
        otfMeanPhase = fftshift(angle(theWaveOTF));
    
        if strcmp(opticsModel, 'WvfHumanMeanOTFmagMeanOTFphase')
           phaseToUse = otfMeanPhase;
        elseif strcmp(opticsModel, 'WvfHumanMeanOTFmagZeroOTFphase')
           phaseToUse = 0*otfMeanPhase;
        end
        
        % Make radially symmetric magnitude
        centerPosition = floor(size(otfMeanMag,1)/2) + 1;
        otfSlice = squeeze(otfMeanMag(centerPosition, :));
        if (waveIndex == 1)
            radiusMatrix = MakeRadiusMat(length(otfSlice),length(otfSlice),centerPosition);
        end
        otfMeanMag = interp1(0:floor(length(otfSlice)/2),otfSlice(centerPosition:end),radiusMatrix,'linear',0);
        customMTFdata.otf(:,:,waveIndex) = ifftshift(otfMeanMag .* exp(1i*phaseToUse));
    end
end

function [otf, translationVector] = otfWithZeroCenteredPSF(otf, psf, translationVector)
    if (isempty(translationVector))
        % Compute center of mass
        centerOfMass = computeCenterOfMass(psf);
   
        % Compute translation vector to bring center of mass at 0,0
        centerPosition = floor(size(psf,1)/2) + 1 * [1 1];
        translationVector = (centerPosition-centerOfMass);
    end
    otf = shiftInFTplane(otf, translationVector);

end

function centerOfMass = computeCenterOfMass(psf0)
    [rc,cc] = ndgrid(1:size(psf0,1),1:size(psf0,2));
    Mt = sum(psf0(:));
    centerOfMassY = sum(psf0(:) .* rc(:)) / Mt;
    centerOfMassX = sum(psf0(:) .* cc(:)) / Mt;
    centerOfMass = [centerOfMassX centerOfMassY];
end

function theWVF = makeWVF(zcoeffs, wavelengthsToCompute, pupilDiameterMM, name)
    theWVF = wvfCreate(...
                'wave',wavelengthsToCompute,...
                'zcoeffs',zcoeffs,...
                'name', name);
    theWVF = wvfSet(theWVF,'calc pupil size',pupilDiameterMM);
    theWVF = wvfComputePSF(theWVF);
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

