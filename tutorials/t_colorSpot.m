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
%   responseParamsGenerate('spatialType','spot','backgroundType','AO','modulationType','AO')
% That function and its subfunctions also documents what the relavant parameters are.
%
% The code illustrated here is encapsulated into function
%   colorSceneCreate.
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
    rParams = responseParamsGenerate('spatialType','spot','backgroundType','AO','modulationType','AO');
end

% Override some defaults to make more sense for our spot application
rParams.oiParams.pupilDiamMm = 7;

%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
theProgram = mfilename;
paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, rParams.colorModulationParams};

%% Make the grayscale spot pattern and have a look
spotPattern = drawSpot(rParams.spatialParams);

% We can see it as a grayscale image
vcNewGraphWin; imagesc(spotPattern); colormap(gray); axis square

% And plot a slice through the center.
%
% This is useful for verifying that the spatial parameters produce the desired
% result in degrees.
figure; hold on;
set(gca,'FontSize',rParams.plotParams.axisFontSize);
xDegs = linspace(-rParams.spatialParams.backgroundSizeDegs/2,rParams.spatialParams.backgroundSizeDegs/2,rParams.spatialParams.col);
plot(xDegs,spotPattern(rParams.spatialParams.row/2,:));
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
nBgWavelengths = length(rParams.backgroundParams.backgroundWavelengthsNm);
bgRadiance = zeros(nWls,1);
for ww = 1:nBgWavelengths
    theWavelength = rParams.backgroundParams.backgroundWavelengthsNm(ww);
    theCornealIrradiance = rParams.backgroundParams.backgroundCornealPowerUW(ww);
    
    % UW is really UW/cm2 because the area of the detector is 1 cm2.  This
    % conversion gives us radiance in UW/[sr-cm2] for the narrowband laser
    % light.
    bgRadianceRaw(ww) = CornIrradianceAndDegrees2ToRadiance(theCornealIrradiance,rParams.spatialParams.backgroundSizeDegs^2);  
    
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
nSpotWavelengths = length(rParams.colorModulationParams.spotWavelengthNm);
spotRadiance = zeros(nWls,1);
for ww = 1:nSpotWavelengths
    theWavelength = rParams.colorModulationParams.spotWavelengthNm(ww);
    theCornealIrradiance = rParams.colorModulationParams.spotCornealPowerUW(ww);
    
    % UW is really UW/cm2 because the area of the detector is 1 cm2.  This
    % conversion gives us radiance in UW/[sr-cm2] for the narrowband laser
    % light.
    spotRadianceRaw(ww) = CornIrradianceAndDegrees2ToRadiance(theCornealIrradiance,rParams.spatialParams.backgroundSizeDegs^2);  
    
    % Convert to Watts/[sr-m2-nm] where we take the wavelength sampling
    % into account so that in the end the calculation of cone responses
    % will come out correctly.
    index = find(theWavelength == wls);
    if (length(index) ~= 1)
        error('Something funky about wls');
    end
    spotRadiance(index) = (10^4)*(10^-6)*spotRadianceRaw(ww)/deltaWl; 
end

%% Produce the isetbio scene

% Create an empty scene to use for the background
sceneBg = sceneCreate('empty');
sceneBg = sceneSet(sceneBg,'wavelength',wls);
sceneBg = sceneSet(sceneBg, 'h fov', rParams.spatialParams.backgroundSizeDegs);

%% Make an image with the background spectral radiance at all locations
radianceEnergyBg = zeros(rParams.spatialParams.row,rParams.spatialParams.col,nWls);
for i = 1:rParams.spatialParams.row
    for j = 1:rParams.spatialParams.col
        radianceEnergyBg(i,j,:) = bgRadiance;
    end
end

%% Convert to quantal units
radiancePhotonsBg = Energy2Quanta(wls,radianceEnergyBg);

%% Put in the photons and the illuminant
%
% This now makes the implied surface reflectance 
% what we started with, as we check a little further
% down.
sceneBg = sceneSet(sceneBg,'photons',radiancePhotonsBg);
%scene = sceneSet(scene,'illuminant energy',theIlluminant);

%% Look at the image contained in our beautiful scene 
%
% This will replace what was in the first figure we had
% and should look the same.
sceneShowImage(sceneBg);

%% Repeat the previous sections for the spot representation
% Creat an empty scene to use for the spot
spotScene = sceneCreate('empty');
spotScene = sceneSet(spotScene,'wavelength',wls);
spotScene = sceneSet(spotScene, 'h fov', rParams.spatialParams.backgroundSizeDegs);

%% Make an image with the background + spot spectral radiance at all locations
radianceEnergySpot = zeros(rParams.spatialParams.row,rParams.spatialParams.col,nWls);
% Background pixels are 1, spot pixels are 2
for i = 1:rParams.spatialParams.row
    for j = 1:rParams.spatialParams.col
        if spotPattern(i,j)== 1; 
            radianceEnergySpot(i,j,:) = bgRadiance;
        elseif spotPattern(i,j) == 2; % Stimulus pixels are flagged "1"
            radianceEnergySpot(i,j,:) = bgRadiance+spotRadiance;
        end
    end
end

%% Convert to quantal units
radiancePhotonsSpot = Energy2Quanta(wls,radianceEnergySpot);

%% Put in the photons and the illuminant
%
% This now makes the implied surface reflectance 
% what we started with, as we check a little further
% down.
spotScene = sceneSet(spotScene,'photons',radiancePhotonsSpot);
%scene = sceneSet(scene,'illuminant energy',theIlluminant);

%% Look at the image contained in our beautiful scene 
%
% This will replace what was in the first figure we had
% and should look the same.
figure;
sceneShowImage(spotScene);

% Look at the scene image.  It is plausible for an L-M grating.  Remember that we
% are looking at the stimuli on a monitor different from the display file
% that we loaded, and thus the RGB values will not produce exactly the
% desired appearance.
vcNewGraphWin; [~,h] = scenePlot(spotScene,'radiance image no grid');
rwObject.write('colorSpotScene',h,paramsList,theProgram,'Type','figure');

%% Create oi
spotOI = oiCreate('wvf human');
spotOI = oiSet(spotOI,'h fov',rParams.spatialParams.backgroundSizeDegs);

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
spotConeMosaic.setSizeToFOV(rParams.spatialParams.backgroundSizeDegs);

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
coneTypes = {'L' 'M' 'S'};
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

