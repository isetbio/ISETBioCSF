function [theCustomOI, Zcoeffs] = oiWithCustomOptics(opticsModel, wavefrontSpatialSamples, calcPupilDiameterMM, umPerDegree, varargin)
    p = inputParser;
    p.addParameter('includeMirrorSubject', false, @islogical);
    p.addParameter('showTranslation', true, @islogical);
    p.addParameter('showMirrorPair', false, @islogical);
    p.parse(varargin{:});
    
    includeMirrorSubject = p.Results.includeMirrorSubject;
    showTranslation = p.Results.showTranslation;
    showMirrorPair = p.Results.showMirrorPair;
    
    % Save rng state. 
    currentState = rng;
    
    % Set rng seed to 1 to ensure we always get the same Zcoeffs
    rng(1);

    if (contains(opticsModel, 'AOoptics'))
        subjectsNum = 1;
        subjectIndices = 1;
        centerPSF = false;
    else     
        switch opticsModel
            case {'WvfHumanMeanOTFmagMeanOTFphase', 'WvfHumanMeanOTFmagZeroOTFphase'}
                subjectsNum = 500;
                subjectIndices = 1:subjectsNum;
                centerPSF = true;

            case {'WvfHuman'}
                subjectsNum = 500;
                subjectIndices = 1;
                centerPSF = true;

            case 'WvfHumanSubject1'
                 subjectsNum = 100;
                 subjectIndices = 1;
                 centerPSF = true;

            case 'WvfHumanSubject2'
                 subjectsNum = 100;
                 subjectIndices = 2; %17;
                 centerPSF = true;

            case 'WvfHumanSubject3'
                 subjectsNum = 100;
                 subjectIndices = 3; %55;
                 centerPSF = true;

            case 'WvfHumanSubject4'
                 subjectsNum = 100;
                 subjectIndices = 4; %69;
                 centerPSF = true;

            case 'WvfHumanSubject5'
                 subjectsNum = 100;
                 subjectIndices = 5; %78;   
                 centerPSF = true;

            otherwise
                error('Unknown optics Model: ''%s''.', opticsModel);
        end
    end
    
    [theCustomOI, Zcoeffs] = generateMeanAcrossSubjectsOTFcustomOI(opticsModel, wavefrontSpatialSamples, calcPupilDiameterMM, subjectsNum, subjectIndices, umPerDegree, centerPSF, includeMirrorSubject, showTranslation, showMirrorPair);
    
    % Restore current rng state
    rng(currentState);
end


function [theCustomOI, Zcoeffs] = generateMeanAcrossSubjectsOTFcustomOI(opticsModel, wavefrontSpatialSamples, calcPupilDiameterMM, subjectsNum, subjectIndices, umPerDegree, centerPSF, includeMirrorSubject, showTranslation, showMirrorPair)

    theOI = oiCreate('wvf human', calcPupilDiameterMM);
    optics = oiGet(theOI,'optics');
    wavelengthsListToCompute = opticsGet(optics,'wave');
   
    % Load the Thibos data set.  If we load 7.5 mm, we can compute for anything possible.
    measPupilDiameterMM = 7.5;
    [Zcoeffs_SampleMean, Zcoeffs_S] = wvfLoadThibosVirtualEyes(measPupilDiameterMM);
    
    Zcoeffs = ieMvnrnd(Zcoeffs_SampleMean, Zcoeffs_S, subjectsNum)'; 
    Zcoeffs = Zcoeffs(:,subjectIndices);
        
    if contains(opticsModel, 'AOoptics')
        Zcoeffs = Zcoeffs * 0;
        measPupilDiameterMM = calcPupilDiameterMM;
    end

    if strcmp(opticsModel, 'WvfHuman')
        availableThibosMeasPupilSizes = [3 4.5 6 7.5];
        idx = find(availableThibosMeasPupilSizes>= calcPupilDiameterMM);
        if (isempty(idx))
            error('There are no Thibos data for pupil size large enough for %2.2f mm', calcPupilDiameterMM);
        end
        measPupilDiameterMM = availableThibosMeasPupilSizes(idx(1));
        [Zcoeffs_SampleMean, Zcoeffs_S] = wvfLoadThibosVirtualEyes(measPupilDiameterMM);
        Zcoeffs = Zcoeffs_SampleMean;
    end
    
    tic
    for subjectIndex = 1:numel(subjectIndices)
        fprintf('Generating wavefront object for ''%s'' optics (subject %d of %d)\n', opticsModel, subjectIndex,numel(subjectIndices));
        
        % Make a WVF object for this subject      
        wvfSubject = makeWVF(wavefrontSpatialSamples, Zcoeffs(:,subjectIndex), wavelengthsListToCompute, ...
            measPupilDiameterMM, calcPupilDiameterMM, umPerDegree, opticsModel);
        
        % Compute an optics structure from the wvf object
        optics = oiGet(wvf2oi(wvfSubject),'optics');
        
        if (contains(opticsModel, 'AOoptics'))
            theCustomOI = oiSet(theOI,'optics', optics);
            continue;
        end
        
        % Proceed with centering of the PSF (for the non AOoptics cases)
        if (centerPSF)
            % Alignment wavelength
            alignmentWavelength = 550;
            [~,alignmentWIndex] = min(abs(wavelengthsListToCompute-alignmentWavelength));
            alignmentWavelength = wavelengthsListToCompute(alignmentWIndex);
    
            % Retrieve OTF and PSF at the alignment wavelength
            alignmentWavelengthOTF = wvfGet(wvfSubject, 'otf', alignmentWavelength);
            alignmentWavelengthPSF = wvfGet(wvfSubject, 'psf', alignmentWavelength);
            
            % Retrieve OTF support in cycles/deg
            xSfCyclesPerRetinalMicron = wvfGet(wvfSubject, 'otf support', 'um', alignmentWavelength);
            xSfCyclesDeg = xSfCyclesPerRetinalMicron * wvfGet(wvfSubject,'um per degree');
            ySfCyclesDeg = xSfCyclesDeg;
            [xSfGridCyclesDegGrid,ySfGridCyclesDegGrid] = meshgrid(xSfCyclesDeg,ySfCyclesDeg);
            
            % Compute translation vector for the PSF at the alignment wavelength
            translationVector = [];
            [~, translationVector, ~, ~, ~] = otfWithZeroCenteredPSF(...
                alignmentWavelengthOTF, alignmentWavelengthPSF, ...
                translationVector, xSfGridCyclesDegGrid,ySfGridCyclesDegGrid, ...
                showTranslation);
        else
            fprintf(2,'Not centering PSF\n');
            translationVector = [0 0];
        end
        
        for wIndex = 1:numel(wavelengthsListToCompute)
            theWavelength = wavelengthsListToCompute(wIndex);
            theWaveOTFbeforeCentering = wvfGet(wvfSubject, 'otf', theWavelength);
            theWavePSFbeforeCentering = wvfGet(wvfSubject, 'psf', theWavelength);
            xSfCyclesPerRetinalMicronAtThisWavelength = wvfGet(wvfSubject, 'otf support', 'um', theWavelength);
            if (centerPSF) && (any(xSfCyclesPerRetinalMicronAtThisWavelength ~= xSfCyclesPerRetinalMicron))
                error('OTF support at %2.0f nm is different than in target wavelength (%2.0fnm)', theWavelength, alignmentWavelength);
            end
            
            if (centerPSF)
                % Obtain an OTF in which the corresponding PSF is
                % translated by the translationVector
                [theWaveOTF, ~, ~, ~,~] = otfWithZeroCenteredPSF(theWaveOTFbeforeCentering, theWavePSFbeforeCentering, translationVector,  xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showTranslation);
            else
                theWaveOTF = theWaveOTFbeforeCentering;
            end
                
            if (subjectIndex == 1) && (wIndex == 1)
                if (contains(opticsModel, 'WvfHumanMean')) && (includeMirrorSubject)
                    centeredOTFs = zeros(4*numel(subjectIndices), size(theWaveOTF,1), size(theWaveOTF,2), numel(wavelengthsListToCompute));
                else
                    centeredOTFs = zeros(numel(subjectIndices), size(theWaveOTF,1), size(theWaveOTF,2), numel(wavelengthsListToCompute));
                end
            end
            
            centeredOTFs(subjectIndex, :,:, wIndex) = theWaveOTF;
            
            if (contains(opticsModel, 'WvfHumanMean')) && (includeMirrorSubject)
                [theOTF90, theOTF180, theOTF270] = ...
                    otfWithMirrorSymmetricPSF(theWaveOTF, xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showMirrorPair);
                centeredOTFs(subjectIndex+numel(subjectIndices), :,:, wIndex) = theOTF90;
                centeredOTFs(subjectIndex+2*numel(subjectIndices), :,:, wIndex) = theOTF180;
                centeredOTFs(subjectIndex+3*numel(subjectIndices), :,:, wIndex) = theOTF270;
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

        radiallySymmetricOTF  = false;
        if (radiallySymmetricOTF)
            tmp = fftshift(otfMeanMag);
            % Make radially symmetric magnitude 
            centerPosition = floor(size(tmp,1)/2) + 1; 
            otfSlice = squeeze(tmp(centerPosition, :)); 
            if (waveIndex == 1) 
                radiusMatrix = MakeRadiusMat(length(otfSlice),length(otfSlice),centerPosition); 
            end 
            tmp2 = interp1(0:floor(length(otfSlice)/2),otfSlice(centerPosition:end),radiusMatrix,'linear',0); 
            otfMeanMag = ifftshift(tmp2);
        
            figure(12345)
            subplot(1,2,1);
            imagesc(xSfCyclesDeg, xSfCyclesDeg, tmp)
            set(gca, 'XLim', [-20 20], 'YLim', [-20 20]);
            axis 'square'

            subplot(1,2,2);
            imagesc(xSfCyclesDeg, xSfCyclesDeg, tmp2)
            set(gca, 'XLim', [-20 20], 'YLim', [-20 20]);
            axis 'square'
            pause
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
                'calc pupil size',calcPupilDiameterMM, ...
                'spatialsamples', wavefrontSpatialSamples, ...
                'zcoeffs', zcoeffs,...
                'name', name);
    
    % Now compute the PSF
    theWVF = wvfComputePSF(theWVF);
end

function [theOTF90, theOTF180, theOTF270] = otfWithMirrorSymmetricPSF(theOTF, xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showMirrorPair)
   theOTF = fftshift(theOTF);
   [~,~, thePSF] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,theOTF);
   thePSF90 = rot90(thePSF,1);
   thePSF180 = rot90(thePSF,2);
   thePSF270 = rot90(thePSF,3);
   [~,~,theOTF90] = PsfToOtf([],[],thePSF90);
   [~,~,theOTF180] = PsfToOtf([],[],thePSF180);
   [~,~,theOTF270] = PsfToOtf([],[],thePSF270);
   
   if (showMirrorPair)
       figure(201); clf;
       subplot(2,2,1);
       imagesc(abs(theOTF))
       subplot(2,2,2);
       imagesc(abs(theOTF180));
       
       subplot(2,2,3);
       imagesc(thePSF)
       subplot(2,2,4);
       imagesc(thePSF180);
       colormap(jet(1024))
       drawnow
   end

   theOTF90 = ifftshift(theOTF90);
   theOTF180 = ifftshift(theOTF180);
   theOTF270 = ifftshift(theOTF270);
end

            
function [centeredOTF,  translationVector, centeredPSF, xGridMinutes,yGridMinutes] = otfWithZeroCenteredPSF(OTF, PSF, translationVector,  xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showTranslation)
    if (isempty(translationVector))
        % Compute center of mass
        %centerOfMass = computeCenterOfMass(PSF);
        centerOfMass = computeCenterOfMassNative(PSF);
        centerPosition = floor(size(PSF,1)/2) + 1 * [1 1];
        
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
