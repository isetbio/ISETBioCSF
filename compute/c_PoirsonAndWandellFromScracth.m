function c_PoirsonAndWandellFromScractch
% c_PoirsonAndWandellFromScratch(varargin)
%
% Compute color detection thresholds to replicate the Poirson & Wandell 1996
                         
    spatialParams = struct(...
        'fieldOfViewDegs', 5.0, ...      In P&W 1996, in the constan cycle condition, this was 10 deg (Section 2.2, p 517)
        'cyclesPerDegree', 2.0,...
             'windowType', 'Gaussian', ...
        'gaussianFWHMDegs', 1.9, ...
        'viewingDistance', 0.75, ...
                   'rows', 256, ...
                   'cols', 256 ...
        );
    
    % In the constant cycle condition, the background was xyY= 0.38, 0.39, 536.2 cd/m2
    % Also they say they placed a uniform field to the screen to increase the contrast resolutionso (page 517, Experimental Aparratus section)
    chromaticParams = struct(...
        'backgroundxyY', [0.38 0.39 536.2], ...
        'testConeContrasts', [0.05 -0.05 0.0] ...
        );
    
    contrast = 0.0;
    backgroundScene = generateGaborDisplayScene(spatialParams, chromaticParams, contrast);
    vcNewGraphWin; [~,h] = scenePlot(backgroundScene,'radiance image no grid');
    pause
    contrast = 1.0;
    modulationScene = generateGaborDisplayScene(spatialParams, chromaticParams, contrast);
    
    vcNewGraphWin; [~,h] = scenePlot(modulationScene,'radiance image no grid');
    pause
    
    
    % Generate the sequence of optical images representing the ramping of the stimulus
    theOIsequence = oiSequenceGenerate(backgroundScene, modulationScene, modulationFunction);
    
    % Generate the cone mosaic with eye movements for theOIsequence
    theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, osNoise, integrationTime, osTimeStep, oiTimeAxis, theOIsequence.length);

    % Compute absorptions and photocurrents
    [absorptionsCountSequence, absorptionsTimeAxis, photoCurrentSequence, photoCurrentTimeAxis] = ...
            theConeMosaic.computeForOISequence(theOIsequence, oiTimeAxis, ...
            'currentFlag', true, ...
            'newNoise', true ...
            );
        
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
    sceneGet(displayScene, 'mean luminance')
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