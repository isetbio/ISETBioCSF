% t_colorGaborScene
%
% Create a scene with a color gabor patch with color directions
% specified as L, M, and S cone contrasts.  The scene will produce
% a Gabor with these contrasts on a specified monitor.
%
% The code illustrated here is encapsulated into function
%   colorGaborSceneCreate.
%
% See also colorGaborSceneCreate, t_colorGaborConeAbsorptionMovie
%
% 7/6/16  dhb  Wrote it.

%% Clear
ieInit; clear; close all;

% Add project toolbox to Matlab path
AddToMatlabPathDynamically(fullfile(fileparts(which(mfilename)),'../toolbox'));

% Some font sizes for plots
labelFontSize = 12;
titleFontSize = 14;
axisFontSize = 8;

%% Define parameters of a gabor pattern
%
% This is for a full contrast vertical gabor centered on the image.
%   fieldOfViewDegs - Field of view in degrees, horizontal direction.
%   cyclesPerDegree - Grating cycles per degree.
%   gaussianFWHMDegs - Full width at half max of spatial Gaussian window.
%   row - Row dimension of scene on monitor
%   col - Col dimension of scene on monitor
%   contrast - Contrast specfied relative to coneContrasts.
%   ang - Angle of grating, in radians
%   ph  - Phase of grating, in radians relative to image center
%   coneContrasts - Color direction of grating in cone contrast space
%   backgroundxYY - Colorimetric specification of background, in CIE xyY (cd/m2)
%   monitorFile - Isetbio display description of monitor on which grating is shown.
%   viewingDistance - Viewing distance of observer from monitor in meters.
gaborParams.fieldOfViewDegs = 4;
gaborParams.cyclesPerDegree = 2;
gaborParams.gaussianFWHMDegs = 1.5;
gaborParams.row = 128;
gaborParams.col = 128;
gaborParams.contrast = 1;
gaborParams.ang = 0;
gaborParams.ph = 0;
gaborParams.coneContrasts = [0.05 -0.05 0]';
gaborParams.backgroundxyY = [0.27 0.30 49.8]';
gaborParams.monitorFile = 'CRT-MODEL';
gaborParams.viewingDistance = 0.75;

% Computed parameters.  These convert numbers to a form used by underlying
% routines.
cyclesPerImage = gaborParams.fieldOfViewDegs*gaborParams.cyclesPerDegree;
gaussianStdDegs = FWHMToStd(gaborParams.gaussianFWHMDegs);
gaussianStdImageFraction = gaussianStdDegs/gaborParams.fieldOfViewDegs;
gaborParams.freq = cyclesPerImage;
gaborParams.GaborFlag = gaussianStdImageFraction;

%% Make the gabor pattern and have a look
%
% We can see it as a grayscale image
gaborPattern = imageHarmonic(gaborParams);
vcNewGraphWin; imagesc(gaborPattern); colormap(gray); axis square

% And plot a slice through the center.
%
% This is useful for verifying that the spatial parameters produce the desired
% result in degrees.  If you generate the Gabor for 0 cpd you can see the Gaussian
% profile and verify that the FWHM is in fact the specified number of
% degrees, and if you make the Gaussian window wide you can count cycles
% and make sure they come out right as well.
figure; hold on;
set(gca,'FontSize',axisFontSize);
xDegs = linspace(-gaborParams.fieldOfViewDegs/2,gaborParams.fieldOfViewDegs/2,gaborParams.col);
plot(xDegs,gaborPattern(gaborParams.row/2,:));
xlabel('Position (degrees)','FontSize',labelFontSize);
ylabel('Image Intensity','FontSize',labelFontSize);

%% Convert Gabor to a color modulation specified in cone space
%
% First, make it a modulation around the mean
% This is easy, because imageHarmoic generates the Gabor as a modulation
% around 1.  Subtracting 1 gives us a modulation in the range -1 to 1.
gaborModulation = gaborPattern-1;

% Convert Gabor to a color modulation specified in cone space
%
% This requires a little colorimetry.
%
% Need to load cone fundamentals and XYZ color matching functions to do the
% needed conversions.  Here we'll use the Stockman-Sharpe 2-degree
% fundamentals and the proposed CIE corresponding XYZ functions.  These
% have the advantage that they are an exact linear transformation away from
% one another.
%
% This is the PTB style data, which I (DHB) know like the back of my hand.  There
% is an Isetbio way to do this too, I'm sure.  The factor of 683 in front
% of the XYZ color matching functions brings the luminance units into cd/m2
% when radiance is in Watts/[m2-sr-nm], which are fairly standard units.
whichXYZ = 'xyzCIEPhys2';
theXYZ = load(['T_' whichXYZ]);
eval(['T_XYZ = 683*theXYZ.T_' whichXYZ ';']);
eval(['S_XYZ = theXYZ.S_' whichXYZ ';']);
clear theXYZ

whichCones = 'cones_ss2';
theCones = load(['T_' whichCones]);
eval(['T_cones = 683*theCones.T_' whichCones ';']);
eval(['S_cones = theCones.S_' whichCones ';']);
clear theCones

% Tranform background into cone excitation coordinates. I always like to
% check with a little plot that I didn't bungle the regression.
M_XYZToCones = ((T_XYZ')\(T_cones'))';
T_conesCheck = M_XYZToCones*T_XYZ;
if (max(abs(T_conesCheck(:)-T_cones(:))) > 1e-3)
    error('Cone fundamentals are not a close linear transform of XYZ CMFs');
end

% Convert background to cone excitations
backgroundConeExcitations = M_XYZToCones*xyYToXYZ(gaborParams.backgroundxyY);

% Convert test cone contrasts to cone excitations
testConeExcitations = (gaborParams.coneContrasts .* backgroundConeExcitations);

% Make the color gabor in LMS excitations
gaborConeExcitationsBg = ones(gaborParams.row,gaborParams.col);
gaborConeExcitations = zeros(gaborParams.row,gaborParams.col,3);
for ii = 1:3
    gaborConeExcitations(:,:,ii) = gaborConeExcitationsBg*backgroundConeExcitations(ii) + ...
        gaborModulation*testConeExcitations(ii);
end

% Check that contrasts come out right.  They will be a little
% less than nominal values becuase it's a Gabor, not a sinusoid.
coneTypes = {'L' 'M' 'S'};
for ii = 1:3
    gaborPlane = gaborConeExcitations(:,:,ii);
    theMax = max(gaborPlane(:)); theMin = min(gaborPlane(:));
    actualConeContrasts(ii) = (theMax-theMin)/(theMax+theMin);
    fprintf('Actual absolute %s cone contrast: %0.3f, nominal: % 0.3f\n', coneTypes{ii}, ...
        actualConeContrasts(ii),abs(gaborParams.coneContrasts(ii)));
end

%% And take a look at the LMS image.  This is just a straight rendering of
% LMS and so won't look the right colors, but we can check that it is
% qualitatively correct.
vcNewGraphWin; imagesc(gaborConeExcitations/max(gaborConeExcitations(:))); axis square

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
display = displayCreate(gaborParams.monitorFile);

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

% Convert the gaborConeExcitations image to RGB
[gaborConeExcitationsCalFormat,m,n] = ImageToCalFormat(gaborConeExcitations);
gaborPrimaryCalFormat = M_ConeExcitationsToPrimary*gaborConeExcitationsCalFormat;
gaborPrimary = CalFormatToImage(gaborPrimaryCalFormat,m,n);

% Check that the image is within the monitor gamut.  If the gabor
% represents an actual stimulus produced with an actual monitor, things
% should be OK if both are represented properly in this routine.
maxPrimary = max(gaborPrimaryCalFormat(:));
minPrimary = min(gaborPrimaryCalFormat(:));
fprintf('Maximum linear RGB (primary) value is %0.2f, minimum %0.2f\n',maxPrimary,minPrimary);
if (maxPrimary > 1 || minPrimary < 0)
    error('RGB primary image is out of gamut.  You need to do something about this.');
end

% Gamma correct the primary values, so we can pop them into an isetbio
% scene in some straightforward manner.
nLevels = size(displayGet(display,'gamma'),1);
gaborRGB = round(ieLUTLinear(gaborPrimary,displayGet(display,'inverse gamma')));

% Make a plot of the gamma correction functions.  These should look
% compressive, the inverse of the monitor gamma function.  In the display
% file CRT-MODEL that we are using, the gamma was given as a power function
% with an exponent of 2, just to model something typical.
vcNewGraphWin; hold on
set(gca,'FontSize',10);
theColors = ['r' 'g' 'b'];
for ii = 1:3
    tempPrimary = gaborPrimary(:,:,ii);
    tempRGB = gaborRGB(:,:,ii);
    plot(tempPrimary(:),tempRGB(:),['o' theColors(ii)],'MarkerFaceColor',theColors(ii));
end
xlim([0 1]);
ylim([0 nLevels]);
axis('square');
xlabel('Linear channel value','FontSize',labelFontSize);
ylabel('Gamma corrected DAC settings','FontSize',labelFontSize);
title('Gamma correction','FontSize',titleFontSize);

% Finally, make the actual isetbio scene
% This combines the image we build and the display properties.
gaborScene = sceneFromFile(gaborRGB,'rgb',[],display);
gaborScene = sceneSet(gaborScene, 'h fov', gaborParams.fieldOfViewDegs);

% Look at the scene image.  It is plausible for an L-M grating.  Remember that we
% are looking at the stimuli on a monitor different from the display file
% that we loaded, and thus the RGB values will not produce exactly the
% desired appearance.
vcNewGraphWin; scenePlot(gaborScene,'radiance image no grid');

%% Create no optics oi
% Create an optical image with no optics.  Here our goal is to verify that
% the cone contrast we compute from actual absorptions matches what we
% specified at the start, so that we are confident we've produced the right
% stimulus scene.  So we don't want to get confused by the optics.
gaborOI = oiCreate('human');
gaborOI = oiSet(gaborOI,'h fov',gaborParams.fieldOfViewDegs);

% Turn off optics for current purpose of checking LMS contrast 
% This involves turning off the off-axis falloff in intensity and
% replacing the OTF with a unity OTF.
optics = oiGet(gaborOI,'optics');
optics = opticsSet(optics,'off axis method','skip');
optics = opticsSet(optics,'OTF',ones(size(opticsGet(optics,'OTF'))));
gaborOI= oiSet(gaborOI,'optics',optics);

% Look at the OI
% Note how different the color appearance is than the scene.  This is
% because the OI incorprates the transmittance of the lens.  Down below we
% will turn that off as a check.
gaborOI = oiCompute(gaborOI,gaborScene);
vcAddAndSelectObject(gaborOI); oiWindow;

%% Verify that removing lens transmittance has expected effect
lens = oiGet(gaborOI,'lens');
lens.density = 0;
gaborOINoLens = oiSet(gaborOI,'lens',lens);
gaborOINoLens = oiCompute(gaborOINoLens,gaborScene);
vcAddAndSelectObject(gaborOINoLens); oiWindow;

%% Create and get noise free sensor using coneMosaic obj
% Create a coneMosaic object here. When setting the fov, if only one value
% is specified, it will automatically make a square cone mosaic.  We don't
% need the whole field of view to check the contrast, which will be
% determined near the center of the Gabor, and making it smaller speeds
% things up.
gaborConeMosaic = coneMosaic;
gaborConeMosaic.setSizeToFOV(gaborParams.fieldOfViewDegs/2);

% There is also an option of whether the cone current should be calculated
% in the compute function. If set to true, it uses an os object inside the
% coneMosaic object. The default is the linearOS.  Here we don't need that.
gaborConeMosaic.noiseFlag = false;
isomerizations = gaborConeMosaic.compute(gaborOI,'currentFlag',false);

%% Take a look at the mosaic responses
gaborConeMosaic.window;

%% Get min max for LMS cone absorptions
% Extract the min and max absorptions in a loop. Since we are
% extracting only L, M, or S absorptions at each iteration, we get a vector
% so one call to max/min will suffice.
%
% These are close enough to the desired values that it seems OK, although
% why they aren't exactly the desired values is a little mysterious.  It is
% possible that the oi/coneMosaic object leads to slightly different cone
% fundamentals, or that the monitor quantization leads to the small
% deviatoins.  Someone energetic could track this down.
conePattern = gaborConeMosaic.pattern;
for ii = 2:4
    maxAbsorption = max(isomerizations(conePattern==ii));
    minAbsorption = min(isomerizations(conePattern==ii));
    
    fprintf('%s cone absorptions\n\tMax: %d \n\tMin: %d\n',coneTypes{ii-1},maxAbsorption,minAbsorption);
    fprintf('\tAbsolute contrast: %04.3f\n',(maxAbsorption-minAbsorption)/(maxAbsorption+minAbsorption));
end
