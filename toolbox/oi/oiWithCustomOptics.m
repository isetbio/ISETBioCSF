function [theCustomOI, Zcoeffs, pupilFunctionData] = oiWithCustomOptics(opticsModel, wavefrontSpatialSamples, calcPupilDiameterMM, umPerDegree, varargin)
    p = inputParser;
    p.addParameter('showTranslation', false, @islogical);
    p.parse(varargin{:});
    
    showTranslation = p.Results.showTranslation;
    
    % Save rng state. 
    currentState = rng;
    
    % Set rng seed to 1 to ensure we always get the same Zcoeffs
    rng(1);

    % Generate default OI
    theOI = oiCreate('wvf human', calcPupilDiameterMM,[],[], umPerDegree);
    optics = oiGet(theOI,'optics');
    
    wavelengthsListToCompute = opticsGet(optics,'wave');
    
    if (strfind(opticsModel, 'AOoptics'))
        fprintf('Generating wavefront object for ''%s'' optics \n', opticsModel);
        [Zcoeffs_SampleMean, ~, ~] = wvfLoadThibosVirtualEyes(7.5);
        Zcoeffs = Zcoeffs_SampleMean * 0;
        measPupilDiameterMM = calcPupilDiameterMM;
        % Generate the WVF    
        [wvfSubject, pupilFunctionData] = makeWVF(wavefrontSpatialSamples, Zcoeffs, wavelengthsListToCompute, ...
            measPupilDiameterMM, calcPupilDiameterMM, umPerDegree, opticsModel);
        % Compute an optics structure from the WVF
        optics = oiGet(wvf2oi(wvfSubject),'optics');
        theCustomOI = oiSet(theOI,'optics', optics);   
        return;
    end
    
    % Load appropriate Thibos dataset
    availableThibosMeasPupilSizes = [3 4.5 6 7.5];
    idx = find(availableThibosMeasPupilSizes>= calcPupilDiameterMM);
    if (isempty(idx))
        error('There are no Thibos data for pupil size large enough for %2.2f mm', calcPupilDiameterMM);
    end
    measPupilDiameterMM = availableThibosMeasPupilSizes(idx(1));
    [Zcoeffs_SampleMean, Zcoeffs_S, subject_coeffs] = wvfLoadThibosVirtualEyes(measPupilDiameterMM);
    ZcoeffSubjects = subject_coeffs.bothEyes;
       
    switch opticsModel
        case {'WvfHumanMeanOTFmagMeanOTFphase', 'WvfHumanMeanOTFmagZeroOTFphase'}
            subjectsNum = 100;
            Zcoeffs = ieMvnrnd(Zcoeffs_SampleMean, Zcoeffs_S, subjectsNum)'; 

        case 'WvfHuman'
            Zcoeffs = Zcoeffs_SampleMean;
            
        case 'ThibosBestPSFSubject3MMPupil'
            Zcoeffs = ZcoeffSubjects(:,98);

        case 'ThibosDefaultSubject3MMPupil'
            Zcoeffs = ZcoeffSubjects(:,132);
                
        case 'ThibosAverageSubject3MMPupil'
            Zcoeffs = ZcoeffSubjects(:,54);
                
        case 'ThibosDefocusedSubject3MMPupil'
            Zcoeffs = ZcoeffSubjects(:,194);
                
        case 'ThibosVeryDefocusedSubject3MMPupil'
            Zcoeffs = ZcoeffSubjects(:,21);
                
        otherwise
            error('Unknown optics Model: ''%s''.', opticsModel);
    end

    
    % Re-align the PSF at 550
    centeringWavelength = 550;
    
    subjectsNum = size(Zcoeffs,2);
    for subjectIndex = 1:subjectsNum
        if (subjectsNum>1)
            fprintf('Generating wavefront object for ''%s'' optics (subject %d of %d)\n', opticsModel, subjectIndex,subjectsNum);
        else
            fprintf('Generating wavefront object for ''%s'' optics \n', opticsModel);
        end
        
        [thePSF, theOTF, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes, pupilFunctionData] = computePSFandOTF(...
            squeeze(Zcoeffs(:,subjectIndex)), ...
            wavelengthsListToCompute, wavefrontSpatialSamples, calcPupilDiameterMM, ...
            centeringWavelength, showTranslation);
        
        if (subjectIndex == 1)
            theMeanOTF = theOTF;
        else
            theMeanOTF = theMeanOTF + theOTF;
        end
    end % subjectIndex
    
    if (subjectsNum>1)
        meanOTF = theMeanOTF/subjectsNum;
        
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
            meanOTF(:,:,waveIndex) = ifftshift(theWaveOTF);
        end
    
        theOTF = meanOTF;
    else
        for waveIndex = 1:numel(wavelengthsListToCompute)
            theWaveOTF = squeeze(theOTF(:,:,waveIndex));
            theOTF(:,:,waveIndex) = ifftshift(theWaveOTF);
        end
    end
    
    % Update optics with new OTF data
    xSfCyclesPerMM = 1000*xSfCyclesDeg / umPerDegree;
    ySfCyclesPerMM = 1000*ySfCyclesDeg / umPerDegree;
    customOptics = opticsSet(optics,'otf data',theOTF);
    customOptics = opticsSet(customOptics, 'otffx',xSfCyclesPerMM);
    customOptics = opticsSet(customOptics,'otffy',ySfCyclesPerMM);
    
    % Update customOI with custom optics
    theCustomOI = oiSet(theOI,'optics', customOptics);
        
    % Restore current rng state
    rng(currentState);
end
