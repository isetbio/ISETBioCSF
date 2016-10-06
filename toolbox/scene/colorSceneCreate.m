function [theScene, gamutScaleFactor] = colorSceneCreate(spatialParams,backgroundParams,colorModulationParams,oiParams,gamutCheckFlag)
% [theScene,gamutScaleFactor] = colorSceneCreate(spatialParams,backgroundParams,colorModulationParams,oiParams,[gamutCheckFlag])
%
% Creates a colored Gabor IBIO scene. The scene will produce a specified
% set of L, M, and S contrasts on a specific monitor.
%
% Inputs:
%   spatialParams    -  A struct that specifies the parameters for the pattern to generate.
%   backgroundParams - A struct that specifies background and display parameters
%   colorModulationParams -  A struct that specifies the parameters for the color modulation.
%   oiParams - A struct that describes the optics.  Here we need (sometimes) it to get
%     the pupil diameter, which we use to back out from corneal irradiance to
%     scene radiance.  Not all threads through this code need this, and it
%     can be empty if it is not needed.
%   gamutCheckFlag  -  If set, the routine returns the scale factor required to bring
%     the cone contrasts into the gamut of the monitor.  This
%     can be useful for setting up simulated experimental
%     conditions.  In this case, the returned scene is
%     empty.  When gamutCheckFlag is false (default), the
%     returned scale factor is 1.
%
% See also t_colorGabor, t_colorSpot, imageHarmonic, drawSpot

%% Optional arg for when we are maximizing contrast
if (nargin < 5 || isempty(gamutCheckFlag))
    gamutCheckFlag = false;
end

%% Extract this field since it's used throughout the function.
%
% All spatial pattern parameter types should define this type
if (~isfield(spatialParams,'fieldOfViewDegs'));
    error('Spatial parameters must have a fieldOfViewDegs field');
end
fieldOfViewDegs = spatialParams.fieldOfViewDegs;

% Background and color modulation parameters must be of compatible type
% Check.  We then assume equality throughout below and only condition on
% colorModulationParams type.
if (~strcmp(colorModulationParams.modulationType,backgroundParams.backgroundType))
    error('Color modulation and background parameters types must be the same');
end

%% Make the spatial pattern
switch(spatialParams.spatialType)
    case 'Gabor'  
        % Check that color modulation/background type is one we understand for Gabors
        switch (colorModulationParams.modulationType)
            case 'monitor'
            case 'AO'
                error('Cannot do a Gabor with AO color modulation/background parameters');
            otherwise
                error('Unknown color modulation/background type specified');
        end
        
        % Make the spatial pattern as a Gabor and convert to a modulation around the mean (mean == 0)
        spatialPattern = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams,colorModulationParams));
        spatialModulation = spatialPattern-1;
        
    case 'spot'
        % Check that color modulation/background type is one we understand
        % for spots
        switch (colorModulationParams.modulationType)
            case 'monitor'
                error('Cannot currently do a spot with monitor-based color modulation/background parameters');
            case 'AO'
            otherwise
                error('Unknown color modulation/background type specified');
        end
        
        % Make the grayscale spot pattern
        spotPattern = drawSpot(spatialParams);
        
    otherwise
        error('Unknown spatial type specified');
end

%% Take background and color modulation into account
switch (colorModulationParams.modulationType)
    case 'monitor'
        % Make sure that the contrast and background vectors are
        % both column vectors.
        coneContrast = colorModulationParams.coneContrasts(:);
        backgroundxyY = backgroundParams.backgroundxyY(:);
        backgroundxyY(3) = backgroundxyY(3)*backgroundParams.lumFactor;
        
        % Convert pattern to a color modulation specified in cone space
        %
        % This requires a little colorimetry.
        %
        % Need to load cone fundamentals and XYZ color matching functions to do the
        % needed conversions.  Here we'll use the Stockman-Sharpe 2-degree
        % fundamentals and the proposed CIE corresponding XYZ functions.  These
        % have the advantage that they are an exact linear transformation away from
        % one another.
        %
        % This is the PTB style data, which I know like the back of my hand.  There
        % is an isetbio way to do this too, I'm sure.  The factor of 683 in front
        % of the XYZ color matching functions brings the luminance units into cd/m2
        % when radiance is in Watts/[m2-sr-nm], which are fairly standard units.
        whichXYZ = 'xyzCIEPhys2';
        theXYZ = load(['T_' whichXYZ]); %#ok<NASGU>
        eval(['T_XYZ = 683*theXYZ.T_' whichXYZ ';']);
        eval(['S_XYZ = theXYZ.S_' whichXYZ ';']);
        clear theXYZ
        
        whichCones = 'cones_ss2';
        theCones = load(['T_' whichCones]); %#ok<NASGU>
        eval(['T_cones = 683*theCones.T_' whichCones ';']);
        eval(['S_cones = theCones.S_' whichCones ';']);
        clear theCones
        
        % Tranform background into cone excitation coordinates. I always like to
        % check with a little plot that I didn't bungle the regression.
        M_XYZToCones = ((T_XYZ')\(T_cones'))'; %#ok<NODEF>
        T_conesCheck = M_XYZToCones*T_XYZ;
        if (max(abs(T_conesCheck(:)-T_cones(:))) > 1e-3)
            error('Cone fundamentals are not a close linear transform of XYZ CMFs');
        end
        
        % Convert background to cone excitations
        backgroundConeExcitations = M_XYZToCones*xyYToXYZ(backgroundxyY);
        
        % Convert test cone contrasts to cone excitations
        testConeExcitationsDir = (coneContrast .* backgroundConeExcitations);
        
        % Make the color pattern in LMS excitations
        patternConeExcitationsBg = ones(spatialParams.row,spatialParams.col);
        patternConeExcitations = zeros(spatialParams.row,spatialParams.col,3);
        for ii = 1:3
            patternConeExcitations(:,:,ii) = patternConeExcitationsBg*backgroundConeExcitations(ii) + ...
                spatialModulation*testConeExcitationsDir(ii);
        end
        
        % Produce an isetbio scene
        %
        % This should represent a monitor image that produces the desired LMS
        % excitations.
        
        % We need a display.  We'll just use the description of a CRT that we have
        % handy.  In doing so, we are assuming that the differences between CRT's
        % used in different threshold experiments do not have a substantial effect
        % on the thresholds.  There will be a little effect because differences in
        % channel spectra will lead to differences in the retinal image because of
        % chromatic aberration, but given the general similarity of monitor channel
        % spectra we expect these differences to be small.  We could check this by
        % doing the calculations with different monitor descriptions.
        display = displayCreate(backgroundParams.monitorFile);
        display = displaySet(display,'spd',backgroundParams.lumFactor*displayGet(display,'spd'));
        
        % Set the viewing distance
        display = displaySet(display,'viewingdistance', spatialParams.viewingDistance);
        
        % Get display channel spectra.  The S vector displayChannelS is PTB format
        % for specifying wavelength sampling: [startWl deltaWl nWlSamples],
        displayChannelWavelengths = displayGet(display,'wave');
        displayChannelS = WlsToS(displayChannelWavelengths);
        displayChannelSpectra = displayGet(display,'spd');
        
        % Spline XYZ and cones to same wavelength sampling as display
        T_conesForDisplay = SplineCmf(S_cones,T_cones,displayChannelWavelengths);
        % T_XYZForDisplay = SplineCmf(S_XYZ,T_XYZ,displayChannelWavelengths);
        
        % Find the matrix that converts between linear channel weights (called
        % "primary" in PTB lingo) and LMS excitations, and its inverse.  Multiplication by
        % the deltaWl is to handle fact that in isetbio radiance is specified in
        % Watts/[m2-sr-nm].
        %
        % Also get matrices for going in and out of XYZ, and compute display max
        % luminance as a sanity check.
        M_PrimaryToConeExcitations = T_conesForDisplay*displayChannelSpectra*displayChannelS(2);
        M_ConeExcitationsToPrimary = inv(M_PrimaryToConeExcitations);
        
        % Find scalar that puts modulation into gamut
        %
        % If we are doing this, we just get the scale factor and return the empty
        % matrix for the scene
        if (gamutCheckFlag)
            backgroundPrimary = M_ConeExcitationsToPrimary*backgroundConeExcitations;
            coneExcitationsDirPrimary = M_ConeExcitationsToPrimary*testConeExcitationsDir;
            gamutScaleFactor = MaximizeGamutContrast(coneExcitationsDirPrimary,backgroundPrimary);
            theScene = [];
            return;
        else
            gamutScaleFactor = 1;
        end
        
        % Convert the cone excitations image to RGB
        [patternConeExcitationsCalFormat,m,n] = ImageToCalFormat(patternConeExcitations);
        patternPrimaryCalFormat = M_ConeExcitationsToPrimary*patternConeExcitationsCalFormat;
        patternPrimary = CalFormatToImage(patternPrimaryCalFormat,m,n);
        
        % Check that the image is within the monitor gamut.  If the pattern
        % represents an actual stimulus produced with an actual monitor, things
        % should be OK if both are represented properly in this routine.
        maxPrimary = max(patternPrimaryCalFormat(:));
        minPrimary = min(patternPrimaryCalFormat(:));
        if (maxPrimary > 1 | minPrimary < 0)
            error('RGB primary image is out of gamut.  You need to do something about this.');
        end
        
        % Gamma correct the primary values, so we can pop them into an isetbio
        % scene in some straightforward manner.  It's important to have a lot of
        % steps in the inverse gamma, so that one doesn't truncate very low
        % contrast scenes.  2^20 seems like a lot.
        patternRGB = round(ieLUTLinear(patternPrimary,displayGet(display,'inverse gamma',2^20)));
        
        % Finally, make the actual isetbio scene
        % This combines the image we build and the display properties.
        theScene = sceneFromFile(patternRGB,'rgb',[],display);
        theScene = sceneSet(theScene, 'h fov', fieldOfViewDegs);
        
    case 'AO'
        % Wavelength sampling
        wls = (colorModulationParams.startWl:colorModulationParams.deltaWl:colorModulationParams.endWl)';
        nWls = length(wls);
        
        % Get equivalent spectral radiance of background and spot increment (full on)
        bgRadiance = AOMonochromaticCornealPowerToRadiance(wls,backgroundParams.backgroundWavelengthsNm,backgroundParams.backgroundCornealPowerUW,pupilDiamMm,spatialParams.backgroundSizeDegs^2);
        spotRadiance = AOMonochromaticCornealPowerToRadiance(wls,colorModulationParams.spotWavelengthsNm,colorModulationParams.sCornealPowerUW,pupilDiamMm,spatialParams.backgroundSizeDegs^2);
        
        % Create an empty scene to use for the spot
        theScene = sceneCreate('empty');
        theScene = sceneSet(theScene,'wavelength',wls);
        theScene = sceneSet(theScene, 'h fov', spatialParams.backgroundSizeDegs);
        
        % Make an image with the background + spot spectral radiance at all locations
        radianceEnergySpot = zeros(spatialParams.row,spatialParams.col,nWls);
        for i = 1:spatialParams.row
            for j = 1:spatialParams.col
                % Background pixels are 1, spot pixels are 2
                if (spotPattern(i,j) == 1)
                    radianceEnergySpot(i,j,:) = bgRadiance;
                elseif (spotPattern(i,j) == 2)
                    radianceEnergySpot(i,j,:) = bgRadiance + colorModulationParams.contrast*spotRadiance;
                end
            end
        end
        radiancePhotonsSpot = Energy2Quanta(wls,radianceEnergySpot);
        
        % Put in the photons 
        theScene = sceneSet(theScene,'photons',radiancePhotonsSpot);
        
        % The way we have set things up, maximum contrast is always 1.
        gamutScaleFactor = 1;
        
    otherwise
        error('Unknown color modulation/background type specified');
end



end

