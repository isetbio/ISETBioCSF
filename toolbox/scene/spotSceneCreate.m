function [spotScene] = spotSceneCreate(spotSpatialParams,backgroundParams,spotSpectralParams)
% [spotScene] = spotSceneCreate(spotSpatialParams,backgroundParams,spotSpectralParams)
% 
% Creates a IBIO scene for a spot such as we can produce in our AO system.
%
% Inputs:
%   spotSpatialParams - A struct that specifies the spatial parameters for
%                       the spot scene.
%   backgroundParams - A struct that specifies the spectral composition of
%                      the background.
%   spotSpectralParams -  A struct that specifies the spectral properties of
%                      the spot
%
% See also t_colorGaborResponseGenerationParams, t_colorGabor, imageHarmonic
%
% 7/7/16 xd   adapted from t_colorGabor
% 7/7/16 npc  added viewing distance param
% 8/11/16 dhb add display factor adjustment
% 9/2/16  dhb increase number of steps in inverse gamma call, to avoid
%             display quantization.

%% Optional arg for when we are maximizing contrast
if (nargin < 4 || isempty(gamutCheckFlag))
    gamutCheckFlag = false;
end

% We also want to make sure that the contrast and background vectors are
% both column vectors.
coneContrast = colorModulationParams.coneContrasts(:);
backgroundxyY = backgroundParams.backgroundxyY(:);
backgroundxyY(3) = backgroundxyY(3)*backgroundParams.lumFactor;

%% Define parameters of a gabor pattern
%
% Extract this field since it's used throughout the function.
fieldOfViewDegs = spatialParams.fieldOfViewDegs;

%% Make the gabor pattern and have a look
%
% We can see it as a grayscale image
gaborPattern = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams,colorModulationParams));

%% Convert the Gabor pattern to a modulation around the mean
%
% This is easy, because imageHarmoic generates the Gabor as a modulation
% around 1.  Subtracting 1 gives us a modulation in the range -1 to 1.
gaborModulation = gaborPattern-1;

%% Convert Gabor to a color modulation specified in cone space
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

% Make the color gabor in LMS excitations
gaborConeExcitationsBg = ones(spatialParams.row,spatialParams.col);
gaborConeExcitations = zeros(spatialParams.row,spatialParams.col,3);
for ii = 1:3
    gaborConeExcitations(:,:,ii) = gaborConeExcitationsBg*backgroundConeExcitations(ii) + ...
        gaborModulation*testConeExcitationsDir(ii);
end

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

%% Find scalar that puts modulation into gamut
% 
% If we are doing this, we just get the scale factor and return the empty
% matrix for the scene
if (gamutCheckFlag)
    backgroundPrimary = M_ConeExcitationsToPrimary*backgroundConeExcitations;
    coneExcitationsDirPrimary = M_ConeExcitationsToPrimary*testConeExcitationsDir;
    gamutScaleFactor = MaximizeGamutContrast(coneExcitationsDirPrimary,backgroundPrimary);
    gaborScene = [];
    return;
else
    gamutScaleFactor = 1;
end

%% Convert the gabor cone excitations image to RGB
[gaborConeExcitationsCalFormat,m,n] = ImageToCalFormat(gaborConeExcitations);
gaborPrimaryCalFormat = M_ConeExcitationsToPrimary*gaborConeExcitationsCalFormat;
gaborPrimary = CalFormatToImage(gaborPrimaryCalFormat,m,n);

%% Check that the image is within the monitor gamut.  If the gabor
% represents an actual stimulus produced with an actual monitor, things
% should be OK if both are represented properly in this routine.
maxPrimary = max(gaborPrimaryCalFormat(:));
minPrimary = min(gaborPrimaryCalFormat(:));
if (maxPrimary > 1 | minPrimary < 0)
    error('RGB primary image is out of gamut.  You need to do something about this.');
end

% Gamma correct the primary values, so we can pop them into an isetbio
% scene in some straightforward manner.  It's important to have a lot of
% steps in the inverse gamma, so that one doesn't truncate very low
% contrast scenes.  2^20 seems like a lot.
gaborRGB = round(ieLUTLinear(gaborPrimary,displayGet(display,'inverse gamma',2^20)));

% Finally, make the actual isetbio scene
% This combines the image we build and the display properties.
gaborScene = sceneFromFile(gaborRGB,'rgb',[],display);
gaborScene = sceneSet(gaborScene, 'h fov', fieldOfViewDegs);

end

