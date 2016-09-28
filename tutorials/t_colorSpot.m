function validationData = t_colorSpot(rParams)
% validationData = t_colorSpot(rParams)
%
% Illustrates the basic steps required to calculate cone isomerizations for
% a monochromatic spot on a background, where the key parameters that will
% vary are the size of the spot, the size of the background, and the
% radiometric properties of the spot and the background.  The reason we
% want to do this is so that we can make predictions for how thresholds
% will vary as we change these properties, particularly the size of the
% spot, for various radiometric, mosaic, and eccentricity choices.
%
% If parameters structure is not passed, the routine will use the defaults
% provided by
%   colorSpotResponseParamsGenerate
% That function and its subfunctions also documents what the relavant parameters are.
%
% The code illustrated here is encapsulated into function
%   colorGaborSceneCreate.
%
% The returned validation structure allows this routine to be called from a
% validation script driven by the UnitTest toolbox.
%
% The tutorial produces output according to a scheme controlled by the
% specified IBIOColorDetect rwObject.
%
% See also:
%
% 9/14/16  dhb, wst  Wrote it.

%% Clear
if (nargin == 0)
    ieInit; close all;
end

%% Fix random number generator so we can validate output exactly
rng(1);

%% Get the parameters we need
%
% t_colorGaborResponseGenerationParams returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
if (nargin < 1 | isempty(rParams))
    rParams = colorSpotResponseParamsGenerate;
end

% Override some defaults to make more sense for our spot application
rParams.oiParams.pupilDiamMm = 7;

%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
theProgram = mfilename;
paramsList = {rParams.spotParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, rParams.colorModulationParams};

%% Make the grayscale spot pattern and have a look
%
% The routine imageHarmonicParamsFromGaborParams massages the
% spotParams/colorModulationParams information into the form needed by
% isetbio's imageHarmonic function.
spotPattern = drawSpot(rParams.spotParams);

% We can see it as a grayscale image
vcNewGraphWin; imagesc(spotPattern); colormap(gray); axis square

% And plot a slice through the center.
%
% This is useful for verifying that the spatial parameters produce the desired
% result in degrees.  If you generate the Gabor for 0 cpd you can see the Gaussian
% profile and verify that the FWHM is in fact the specified number of
% degrees, and if you make the Gaussian window wide you can count cycles
% and make sure they come out right as well.
figure; hold on;
set(gca,'FontSize',rParams.plotParams.axisFontSize);
xDegs = linspace(-rParams.spotParams.backgroundSizeDegs/2,rParams.spotParams.backgroundSizeDegs/2,rParams.spotParams.col);
plot(xDegs,spotPattern(rParams.spotParams.row/2,:));
xlabel('Position (degrees)','FontSize',rParams.plotParams.labelFontSize);
ylabel('Image Intensity','FontSize',rParams.plotParams.labelFontSize);

%% Take the measurements made for the AO stimulus and compute the radiance we need

% Wavelengths sampling
startWl = 550;
endWl = 830;
deltaWl = 10;
wls = (startWl:deltaWl:endWl)';
nWls = length(wls);

% Background
nBgWavelenths = length(rParams.spotParams.backgroundWavelengthsNm);
bgRadiance = zeros(nWls,1);
for ww = 1:nBgWavelenths
    theWavelength = rParams.spotParams.backgroundWavelengthsNm(ww);
    theCornealIrradiance = rParams.spotParams.backgroundCornealIrradianceUW(ww);
    
    % UW is really UW/cm2 because the area of the detector is 1 cm2.  This
    % conversion gives us radiance in UW/[sr-cm2] for the narrowband laser
    % light.
    bgRadianceRaw(ww) = CornIrradianceAndDegrees2ToRadiance(theCornealIrradiance,rParams.spotParams.backgroundSizeDegs^2);  
    
    % Convert to Watts/[sr-m2-nm] where we take the wavelength sampling
    % into account so that in the end the calculation of cone responses
    % will come out correctly.
    index = find(theWavelength == wls);
    if (length(index) ~= 1)
        error('Something funky about wls');
    end
    bgRadiance(index) = (10^4)*(10^-6)*bgRadianceRaw(ww)/deltaWl; 
end


% Spot

%% Produce the isetbio scene

% Create an empty scene
scene = sceneCreate('empty');
scene = sceneSet(scene,'wavelength',wls);

%% Make an image with the background spectral radiance at all locations
radianceEnergy = zeros(rParams.spotParams.row,rParams.spotParams.col,nWls);
for i = 1:rParams.spotParams.row
    for j = 1:rParams.spotParams.col
        radianceEnergy(i,j,:) = bgRadiance;
    end
end

%% Convert to quantal units
radiancePhotons = Energy2Quanta(wls,radianceEnergy);

%% Put in the photons and the illuminant
%
% This now makes the implied surface reflectance 
% what we started with, as we check a little further
% down.
scene = sceneSet(scene,'photons',radiancePhotons);
%scene = sceneSet(scene,'illuminant energy',theIlluminant);

%% Look at the image contained in our beautiful scene 
%
% This will replace what was in the first figure we had
% and should look the same.
sceneShowImage(scene);


%% Produce an isetbio scene
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
display = displayCreate(rParams.backgroundParams.monitorFile);
display = displaySet(display,'viewingdistance', rParams.spotParams.viewingDistance);

% Get display channel spectra.  The S vector displayChannelS is PTB format
% for specifying wavelength sampling: [startWl deltaWl nWlSamples],
displayChannelWavelengths = displayGet(display,'wave');
displayChannelS = WlsToS(displayChannelWavelengths);
displayChannelSpectra = displayGet(display,'spd');

% Spline XYZ and cones to same wavelength sampling as display
T_conesForDisplay = SplineCmf(S_cones,T_cones,displayChannelWavelengths);
T_XYZForDisplay = SplineCmf(S_XYZ,T_XYZ,displayChannelWavelengths);

% Find the matrix that converts between linear channel weights (called
% "primary" in PTB lingo) and LMS excitations, and its inverse.  Multiplication by
% the deltaWl is to handle fact that in isetbio radiance is specified in
% Watts/[m2-sr-nm].
%
% Also get matrices for going in and out of XYZ, and compute display max
% luminance as a sanity check.
M_PrimaryToConeExcitations = T_conesForDisplay*displayChannelSpectra*displayChannelS(2);
M_ConeExcitationsToPrimary = inv(M_PrimaryToConeExcitations);

M_PrimaryToXYZ = T_XYZForDisplay*displayChannelSpectra*displayChannelS(2);
M_XYZToPrimary = inv(M_PrimaryToXYZ);
displayMaxXYZ = M_PrimaryToXYZ*[1 1 1]';
fprintf('Max luminace of the display is %0.1f cd/m2\n',displayMaxXYZ(2));

% Convert the spotConeExcitations image to RGB
[spotConeExcitationsCalFormat,m,n] = ImageToCalFormat(spotConeExcitations);
spotPrimaryCalFormat = M_ConeExcitationsToPrimary*spotConeExcitationsCalFormat;
spotPrimary = CalFormatToImage(spotPrimaryCalFormat,m,n);

% Check that the image is within the monitor gamut.  If the spot
% represents an actual stimulus produced with an actual monitor, things
% should be OK if both are represented properly in this routine.
maxPrimary = max(spotPrimaryCalFormat(:));
minPrimary = min(spotPrimaryCalFormat(:));
fprintf('Maximum linear RGB (primary) value is %0.2f, minimum %0.2f\n',maxPrimary,minPrimary);
if (maxPrimary > 1 || minPrimary < 0)
    error('RGB primary image is out of gamut.  You need to do something about this.');
end

% Gamma correct the primary values, so we can pop them into an isetbio
% scene in some straightforward manner.
nLevels = size(displayGet(display,'gamma'),1);
spotRGB = round(ieLUTLinear(spotPrimary,displayGet(display,'inverse gamma')));

% Make a plot of the gamma correction functions.  These should look
% compressive, the inverse of the monitor gamma function.  In the display
% file CRT-MODEL that we are using, the gamma was given as a power function
% with an exponent of 2, just to model something typical.
vcNewGraphWin; hold on
set(gca,'FontSize',10);
theColors = ['r' 'g' 'b'];
for ii = 1:3
    tempPrimary = spotPrimary(:,:,ii);
    tempRGB = spotRGB(:,:,ii);
    plot(tempPrimary(:),tempRGB(:),['o' theColors(ii)],'MarkerFaceColor',theColors(ii));
end
xlim([0 1]);
ylim([0 nLevels]);
axis('square');
xlabel('Linear channel value','FontSize',rParams.plotParams.labelFontSize);
ylabel('Gamma corrected DAC settings','FontSize',rParams.plotParams.labelFontSize);
title('Gamma correction','FontSize',rParams.plotParams.titleFontSize);

% Finally, make the actual isetbio scene
% This combines the image we build and the display properties.
spotScene = sceneFromFile(spotRGB,'rgb',[],display);
spotScene = sceneSet(spotScene, 'h fov', rParams.spotParams.backgroundSizeDegs);

% Look at the scene image.  It is plausible for an L-M grating.  Remember that we
% are looking at the stimuli on a monitor different from the display file
% that we loaded, and thus the RGB values will not produce exactly the
% desired appearance.
vcNewGraphWin; [~,h] = scenePlot(spotScene,'radiance image no grid');
rwObject.write('colorSpotScene',h,paramsList,theProgram,'Type','figure');

%% Create oi
spotOI = oiCreate('wvf human');
spotOI = oiSet(spotOI,'h fov',rParams.spotParams.backgroundSizeDegs);

% Set pupil diameter
focalLength = oiGet(spotOI,'distance');
desiredFNumber = focalLength/(rParams.oiParams.pupilDiamMm/1000);
spotOI  = oiSet(spotOI ,'optics fnumber',desiredFNumber);
pupilDiamMmCheck = 1000*oiGet(spotOI,'optics aperture diameter');
if (max(abs(pupilDiamMmCheck - rParams.oiParams.pupilDiamMm)) > 1e-8)
    error('Failed to set pupil diameter as expected');
end

%% Compute blurred optical image
%
% Turn of default off axis intensity falloff calculation first.
optics = oiGet(spotOI,'optics');
optics = opticsSet(optics,'off axis method','skip');
spotOI = oiSet(spotOI,'optics',optics);
spotOIBlur = oiCompute(spotOI,spotScene);

% Note how different the color appearance is than the scene.  This is
% because the OI incorprates the transmittance of the lens.  Down below we
% will turn that off as a check.
vcNewGraphWin; [~,h] = oiPlot(spotOIBlur,'irradiance image no grid');
rwObject.write('colorGaborOpticalImageBlur',h,paramsList,theProgram,'Type','figure');
clearvars('spotOIBlur');

%% Turn off optics for current purpose of checking LMS contrast 
% This involves replacing the OTF with a unity OTF, and recompute
optics = opticsSet(optics,'OTF',ones(size(opticsGet(optics,'OTF'))));
spotOI = oiSet(spotOI,'optics',optics);
spotOI = oiCompute(spotOI,spotScene);

% Look at the OI
vcNewGraphWin; [~,h] = oiPlot(spotOI,'irradiance image no grid');
rwObject.write('colorGaborOpticalImageNoBlur',h,paramsList,theProgram,'Type','figure');

% Just for fun, put OI into isetbio's interactive window
vcAddAndSelectObject(spotOI); oiWindow;

%% Verify that removing lens transmittance has expected effect
lens = oiGet(spotOI,'lens');
lens.density = 0;
spotOINoLens = oiSet(spotOI,'lens',lens);
spotOINoLens = oiCompute(spotOINoLens,spotScene);
vcAddAndSelectObject(spotOINoLens); oiWindow;

%% Create and get noise free sensor using coneMosaic obj
% Create a coneMosaic object here. When setting the fov, if only one value
% is specified, it will automatically make a square cone mosaic.
spotConeMosaic = coneMosaic;
spotConeMosaic.setSizeToFOV(rParams.spotParams.backgroundSizeDegs);

% There is also an option of whether the cone current should be calculated
% in the compute function. If set to true, it uses an os object inside the
% coneMosaic object. The default is the linearOS.  Here we don't need that.
spotConeMosaic.noiseFlag = false;
isomerizations = spotConeMosaic.compute(spotOI,'currentFlag',false);

%% Take a look at the mosaic responses in the window
spotConeMosaic.window;

% And must make a plot in a figure
vcNewGraphWin; [~,h] = spotConeMosaic.plot('cone mosaic');
rwObject.write('colorGaborMosaic',h,paramsList,theProgram,'Type','figure');
vcNewGraphWin; [~,h] = spotConeMosaic.plot('mean absorptions');
rwObject.write('colorGaborIsomerizations',h,paramsList,theProgram,'Type','figure');

%% Get min max for LMS cone isomerizations
% Extract the min and max absorptions in a loop. Since we are
% extracting only L, M, or S absorptions at each iteration, we get a vector
% so one call to max/min will suffice.
%
% These are close enough to the desired values that it seems OK, although
% why they aren't exactly the desired values is a little mysterious.  It is
% possible that the oi/coneMosaic object leads to slightly different cone
% fundamentals, or that the monitor quantization leads to the small
% deviatoins.  Someone energetic could track this down.
conePattern = spotConeMosaic.pattern;
for ii = 2:4
    maxIsomerizations(ii) = max(isomerizations(conePattern==ii));
    minIsomerizations(ii) = min(isomerizations(conePattern==ii));
    contrasts(ii) = ...
        (maxIsomerizations(ii)-minIsomerizations(ii))/(maxIsomerizations(ii)+minIsomerizations(ii));
    
    fprintf('%s cone isomerzations\n\tMax: %d \n\tMin: %d\n',coneTypes{ii-1},maxIsomerizations(ii),minIsomerizations(ii));
    fprintf('\tAbsolute contrast: %04.3f\n',contrasts(ii));
end

%% Send back some validation data if requested
if (nargout > 0)
    validationData.maxIsomerizations = maxIsomerizations;
    validationData.minIsomerizations = minIsomerizations;
    validationData.contrasts = contrasts;
end

