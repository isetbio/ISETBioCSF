function [theOI, theCustomOI, Zcoeffs, populationOTFdata] = oiWithCustomOptics(pupilDiameterMM, opticsModel)
    % 
    theOI = oiCreate('wvf human', pupilDiameterMM);
    theCustomOI = theOI;

    % Save rng state. 
    currentState = rng;
    
    % Set rng seed to 1 to ensure we always get the same Zcoeffs
    rng(1);
    
    % Load the 4.5 mm pupil diameter Thibos data set
    measPupilDiameterMM = 4.5;
    [Zcoeffs_SampleMean, Zcoeffs_S] = wvfLoadThibosVirtualEyes(measPupilDiameterMM);

    subjectsNum  = 100;
    % Generate Z-coeffs for subjectsNum
    Zcoeffs = ieMvnrnd(Zcoeffs_SampleMean, Zcoeffs_S, subjectsNum)'; 

    % Restore current rng state
    rng(currentState);
    
    % Pull out the optics structure
    optics = oiGet(theCustomOI,'optics');
    wavelengthsListToCompute = opticsGet(optics,'wave');
    % Compute 
    [customOTFdata, Zcoeffs, populationOTFdata] = computeCustomOTF(Zcoeffs, pupilDiameterMM, wavelengthsListToCompute, opticsModel);
    
    % Update optics with new OTF data
    optics = opticsSet(optics,'otf data',customOTFdata.otf);

    % Stick optics into oi
    theCustomOI = oiSet(theCustomOI,'optics',optics);
end

function [customOTFdata, Zcoeffs, populationOTFdata] = computeCustomOTF(Zcoeffs, pupilDiameterMM, wavelengths, opticsModel)

    if ismember(opticsModel, {'WvfHumanMeanOTFmagMeanOTFphase', 'WvfHumanMeanOTFmagZeroOTFphase'})
        [customOTFdata, Zcoeffs,  populationOTFdata] = computeMeanSubjectOTF(Zcoeffs, pupilDiameterMM, wavelengths, opticsModel);
    else
        [customOTFdata, Zcoeffs] = computeSingleSubjectOTF(Zcoeffs, pupilDiameterMM, wavelengths, opticsModel);
        populationOTFdata = [];
    end
end

function [customOTFdata, Zcoeffs] = computeSingleSubjectOTF(Zcoeffs, pupilDiameterMM, wavelengths, opticsModel)
     % Make WVF for this subject   
     switch opticsModel
         case 'WvfHumanSubject1'
             subjectIndex = 1;
         case 'WvfHumanSubject2'
             subjectIndex = 17;
         case 'WvfHumanSubject3'
             subjectIndex = 55;
         case 'WvfHumanSubject4'
             subjectIndex = 69;
         case 'WvfHumanSubject5'
            subjectIndex = 78;
         otherwise
             error('Unknown optics model: ''%s''.', opticsModel);
     end
     
     % Make WVF for this subject  
     Zcoeffs = Zcoeffs(:,subjectIndex);
     wvfSubject = makeWVF(Zcoeffs, wavelengths, pupilDiameterMM, sprintf('subject-%d', subjectIndex));
     optics = oiGet(wvf2oi(wvfSubject),'optics'); 
     
     % Get OTF
     otf = opticsGet(optics, 'otf');
        
     % Get PSF at reference wavelength
     referencePSF = wvfGet(wvfSubject, 'psf', 550);
     
     % Center PSFs
     for waveIndex = 1:numel(wavelengths)
        % Get OTF at  wavelength
        theWaveOTF = squeeze(otf(:,:,waveIndex));
            
        if (waveIndex == 1)
           [theWaveOTF, translationVector] = otfWithZeroCenteredPSF(theWaveOTF, referencePSF, []);
        else
           theWaveOTF = otfWithZeroCenteredPSF(theWaveOTF, referencePSF, translationVector);
        end
        otf(:,:,waveIndex) = theWaveOTF;
     end % waveIndex
     
     customOTFdata.xSfGridCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
     customOTFdata.ySfGridCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
     customOTFdata.otf = otf;
end

function [customOTFdata, Zcoeffs, otfSubjects] = computeMeanSubjectOTF(Zcoeffs, pupilDiameterMM, wavelengths,  opticsModel)
    
    subjectsNum = size(Zcoeffs,2);
    for subjectIndex = 1:subjectsNum
        fprintf('Computing wavefront data for subject %d/%d\n', subjectIndex,subjectsNum);
        % Make WVF for this subject      
        wvfSubject = makeWVF(Zcoeffs(:,subjectIndex), wavelengths, pupilDiameterMM, sprintf('subject-%d', subjectIndex));
        optics = oiGet(wvf2oi(wvfSubject),'optics');
        
        % Get OTF
        otf = opticsGet(optics, 'otf');
        
        % Get PSF at reference wavelength
        referencePSF = wvfGet(wvfSubject, 'psf', 550);

        % Center PSFs
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
                customOTFdata.xSfGridCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
                customOTFdata.ySfGridCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
                customOTFdata.otf = 0*otf;
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
        customOTFdata.otf(:,:,waveIndex) = ifftshift(otfMeanMag .* exp(1i*phaseToUse));
    end
end

function [OTF, translationVector] = otfWithZeroCenteredPSF(OTF, PSF, translationVector)
    if (isempty(translationVector))
        % Compute center of mass
        centerOfMass = computeCenterOfMass(PSF);
   
        % Compute translation vector to bring center of mass at 0,0
        centerPosition = floor(size(PSF,1)/2) + 1 * [1 1];
        translationVector = (centerPosition-centerOfMass);
    end
    OTF = shiftInFTplane(OTF, translationVector);
end

function centerOfMass = computeCenterOfMass(PSF)
    [rc,cc] = ndgrid(1:size(PSF,1),1:size(PSF,2));
    Mt = sum(PSF(:));
    centerOfMassY = sum(PSF(:) .* rc(:)) / Mt;
    centerOfMassX = sum(PSF(:) .* cc(:)) / Mt;
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