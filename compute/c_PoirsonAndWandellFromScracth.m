function c_PoirsonAndWandellFromScractch
% c_PoirsonAndWandellFromScratch(varargin)
%
% Compute color detection thresholds to replicate the Poirson & Wandell 1996
         
    close all
    spatialParams = struct(...
        'fieldOfViewDegs', 5.0, ...      In P&W 1996, in the constan cycle condition, this was 10 deg (Section 2.2, p 517)
        'cyclesPerDegree', 2.0,...
      'spatialPhaseInDeg', 90, ...
             'windowType', 'Gaussian', ...
        'gaussianFWHMDegs', 1.9, ...
        'viewingDistance', 0.75, ...
                    'row', 256, ...
                    'col', 256 ...
        );
    
    CRTrefreshRate = 87;
    temporalParams = struct(...
        'stimSamplingInterval', 1.0/CRTrefreshRate, ...
                 'envelopeTau', 165/1000, ...
               'envelopeWidth', 4.0/CRTrefreshRate, ...
                'envelopePeak', 0/1000 ...
        );
    
    mosaicParams = struct(...
                'size', 0.1*spatialParams.fieldOfViewDegs,...         % nan for 1L, 1M, and 1S-cone only
     'integrationTime', 50/1000, ...        % we will vary this one
         'photonNoise', true, ...           % add Poisson noise
          'osTimeStep', 1/1000, ...         % 1 millisecond
             'osNoise', false ...           % outer-segment noise
       );
    
    
    % In the constant cycle condition, the background was xyY= 0.38, 0.39, 536.2 cd/m2
    % Also they say they placed a uniform field to the screen to increase the contrast resolutionso (page 517, Experimental Aparratus section)
    xyY =  [0.38 0.39 536.2];
    backgroundChromaticParams = struct(...
        'backgroundxyY', xyY, ...
        'testConeContrasts', [0.0 0.0 0.0] ...
        );
    
    stimulusChromaticParams = struct(...
        'backgroundxyY', xyY, ...
        'testConeContrasts', [0.5 0.5 0.5] ...
        );
    
    % Define how the stimulus is modulated in time
    stimTimeAxis = -temporalParams.envelopeWidth:temporalParams.stimSamplingInterval:temporalParams.envelopeWidth;
    % monophasic modulation function
    stimulusModulationFunction = exp(-0.5*((stimTimeAxis-temporalParams.envelopePeak)/temporalParams.envelopeTau).^2);
    
    % Generate background scene
    backgroundScene = generateGaborDisplayScene(spatialParams, backgroundChromaticParams, 0.0);
    
    % Generate modulated scene with testContrast
    testContrast = 1.0;
    modulatedScene = generateGaborDisplayScene(spatialParams, stimulusChromaticParams, testContrast);
    
    % Generate custom optics
    oiParams = struct(...
        'fieldOfViewDegs', 5.0 ...
        );
    theOI = oiGenerate(oiParams);
    
    % Compute background OI
    oiBackground = theOI;
    oiBackground = oiCompute(oiBackground, backgroundScene);
    
    % Compute modulated OI
    oiModulated = theOI;
    oiModulated = oiCompute(oiModulated, modulatedScene);
    
    % Generate the sequence of optical images representing the ramping of the stimulus
    theOIsequence = oiSequence(oiBackground, oiModulated, stimTimeAxis, stimulusModulationFunction, 'composition', 'blend');
    theOIsequence.visualize();
    
    % Generate the cone mosaic with eye movements for theOIsequence
    theConeMosaic = coneMosaicGenerate(mosaicParams, theOIsequence);

    
    
    % Compute absorptions and photocurrents
    [absorptionsCountSequence, absorptionsTimeAxis, photoCurrentSignals, photoCurrentTimeAxis] = ...
            theConeMosaic.computeForOISequence(theOIsequence, ...
            'currentFlag', true, ...
            'newNoise', true ...
            );

    % Visualize oisequence with eye movement sequence
    theOIsequence.visualizeWithEyeMovementSequence(absorptionsTimeAxis);
    
    % Plot responses
    iL = find(theConeMosaic.pattern == 2);
    iM = find(theConeMosaic.pattern == 3);
    iS = find(theConeMosaic.pattern == 4);
    
    plotResponseTimeSeries(absorptionsTimeAxis, absorptionsCountSequence, iL, iM, iS);
    plotResponseTimeSeries(photoCurrentTimeAxis, photoCurrentSignals, iL, iM, iS);
end

function plotResponseTimeSeries(timeAxis, signals, iL, iM, iS)

    % Reshape
    [signals, ~, ~] = RGB2XWFormat(signals);
    
    figure(); clf;
    plot(timeAxis*1000, squeeze(signals(iL,:)), 'r.');
    hold on;
    plot(timeAxis*1000, squeeze(signals(iM,:)), 'g.');
    plot(timeAxis*1000, squeeze(signals(iS,:)), 'b.');
    xlabel('time (ms)');
    title(sprintf('%d L-cones, %d M-cones, %d S-cones', numel(iL), numel(iM), numel(iS)));
end

function theConeMosaic = coneMosaicGenerate(mosaicParams, theOIsequence)
    % Default human mosaic
    theConeMosaic = coneMosaic;
    
    % Adjust size
    if isnan(mosaicParams.size)
        % Generate a human cone mosaic with 1L, 1M and 1S cone
        theConeMosaic.rows = 1;
        theConeMosaic.cols = 3;
        theConeMosaic.pattern = [2 3 4];
    else
        theConeMosaic.setSizeToFOV(mosaicParams.size);
    end
    
    % Set the noise
    theConeMosaic.noiseFlag = mosaicParams.photonNoise;

    % Set the integrationTime
    theConeMosaic.integrationTime = mosaicParams.integrationTime;
    
    % Generate the outer-segment object to be used by the coneMosaic
    theOuterSegment = osLinear();
    theOuterSegment.noiseFlag = mosaicParams.osNoise;
    
    % Set a custom timeStep, for @osLinear we do not need the default 0.1 msec
    theOuterSegment.timeStep = mosaicParams.osTimeStep;

    % Couple the outersegment object to the cone mosaic object
    theConeMosaic.os = theOuterSegment;

    % Generate eye movement sequence for all oi's
    stimulusSamplingInterval = theOIsequence.oiTimeAxis(2)-theOIsequence.oiTimeAxis(1);
    eyeMovementsNumPerOpticalImage = stimulusSamplingInterval/theConeMosaic.integrationTime;
    eyeMovementsNum = round(eyeMovementsNumPerOpticalImage*theOIsequence.length);
    
    if (eyeMovementsNum < 1)
        error('Less than 1 eye movement!!! \nStimulus sampling interval:%g Cone mosaic integration time: %g\n', stimulusSamplingInterval, theConeMosaic.integrationTime);
    else 
        fprintf('Optical image sequence contains %2.0f eye movements (%2.2f eye movements/oi)\n', eyeMovementsNum, eyeMovementsNumPerOpticalImage);
        theConeMosaic.emGenSequence(eyeMovementsNum);
    end
end


function theOI = oiGenerate(oiParams)
    % Generate optics
    if ((isfield(oiParams, 'noOptics')) && (oiParams.noOptics == true))
        theOI = oiCreate('diffraction limited');
        optics = oiGet(theOI,'optics');           
        optics = opticsSet(optics,'fnumber',0);
        optics = opticsSet(optics, 'off axis method', 'skip');
        theOI = oiSet(theOI,'optics', optics);
    else
        theOI = oiCreate('wvf human');
    end
    
    % Remove Lens
    if ((isfield(oiParams, 'noLens')) && (oiParams.noLens == true))
        lens = oiGet(theOI,'lens');
        lens.density = 0;
        theOI = oiSet(theOI,'lens',lens);
    end
    
    % Set custom FOV
    if (isfield(oiParams, 'fieldOfViewDegs'))
        theOI = oiSet(theOI,'h fov',oiParams.fieldOfViewDegs);
    end
    
    % Set custom pupil diameter
    if (isfield(oiParams, 'pupilDiamMM'))
        focalLength = oiGet(theOI,'distance');
        desiredFNumber = focalLength/(oiParams.pupilDiamMm/1000);
        theOI  = oiSet(theOI ,'optics fnumber',desiredFNumber);
        pupilDiamMmCheck = 1000*oiGet(theOI,'optics aperture diameter');
        if (max(abs(pupilDiamMmCheck - oiParams.pupilDiamMm)) > 1e-8)
            error('Failed to set pupil diameter as expected');
        end
    end
    
end


function displayScene = generateGaborDisplayScene(spatialParams, colorParams, contrast)

    % Genereate the rendering display
    display = displayCreate('CRT-MODEL');
    display = displaySet(display,'viewingdistance', 1.0);
    
    % Adjust display's SPDs so as to be able to generate the desired luminance
    peakLuminanceBeforeAdjustment = displayGet(display, 'peak luminance');
    spdMultiplierToGetDesiredLum = colorParams.backgroundxyY(3)/(0.5*peakLuminanceBeforeAdjustment);
    display = displaySet(display,'spd',spdMultiplierToGetDesiredLum*displayGet(display,'spd'));
    
    % Get the SPDs and their spectral sampling
    displayWls = displayGet(display,'wave');
    displaySpectralSampling = WlsToS(displayWls);
    displaySPD = displayGet(display,'spd');
    
    % Load cone fundamentals and XYZtoConeExcitationsMatrix
    [XYZtoConeExcitationsMatrix, coneFundamentals, coneSpectralSampling] = XYZToCones();
    
    % Resample cone fundamentals to the spectral sampling of the display
    coneFundamentals = SplineCmf(coneSpectralSampling, coneFundamentals, displaySpectralSampling);
    
    % Compute coneExcitationsToDisplayPrimaryMatrix
    displayPrimaryToConeExcitations = coneFundamentals * displaySPD * displaySpectralSampling(2);
    coneExcitationsToDisplayPrimaryMatrix = inv(displayPrimaryToConeExcitations);
    
    % Convert the background xyY to cone excitations
    backgroundConeExcitations = XYZtoConeExcitationsMatrix * xyYToXYZ(colorParams.backgroundxyY');
    
    % Generate a Gabor spatial modulation pattern
    gaborModulationPattern = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, contrast))-1;
    
    % Compute cone excitationImage using the gaborModulationPattern, the
    % backgroundConeExcitations and the testConeContrasts
    for ii = 1:3
        coneExcitationImage(:,:,ii) = backgroundConeExcitations(ii) * (1 + gaborModulationPattern * colorParams.testConeContrasts(ii));
    end

    % Compute the RGB primary values
    [coneExcitations,m,n] = ImageToCalFormat(coneExcitationImage);
    displayRGBPrimaries = coneExcitationsToDisplayPrimaryMatrix * coneExcitations;
    displayRGBPrimaryImage = CalFormatToImage(displayRGBPrimaries,m,n);
    
    % Check for out-of-gamut
    if (any(displayRGBPrimaryImage(:) > 1))
        fprintf(2,'Image above gamut\n');
    end
    if (any(displayRGBPrimaryImage(:) < 0))
        fprintf(2,'Image below gamut\n');
    end
    
    % Gamma correct the primary values
    displayRGBsettingsImage = round(ieLUTLinear(displayRGBPrimaryImage,displayGet(display,'inverse gamma')));
    
    % Make a scene from the gaborRGBsettings values
    displayScene = sceneFromFile(displayRGBsettingsImage,'rgb',[],display);
    displayScene = sceneSet(displayScene, 'h fov', spatialParams.fieldOfViewDegs);
    displayScene = sceneSet(displayScene, 'distance', spatialParams.viewingDistance);
end

function [M, coneFundamentals, coneSpectralSampling] = XYZToCones
    % Here we'll use the Stockman-Sharpe 2-degree fundamentals and the proposed CIE corresponding XYZ functions
    theCones = load('T_cones_ss2');
    theXYZ = load('T_xyzCIEPhys2');
    
    XYZcolorMatchingFunctions = 683 * theXYZ.T_xyzCIEPhys2;
    coneFundamentals = 683 * theCones.T_cones_ss2;
    coneSpectralSampling = theCones.S_cones_ss2;
    M = ((XYZcolorMatchingFunctions')\(coneFundamentals'))';
end