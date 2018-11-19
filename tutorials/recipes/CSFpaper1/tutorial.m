%% Tutorial listed in the CSF paper
%% NPC, ISETBIO Team
%%
%% Generate a display for presenting stimuli and place it at a viewing distance of 57 cm
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 0.57);

%% Specify a Gabor stimulus 
stimParams = struct(...
    'spatialFrequencyCyclesPerDeg', 10, ... % 10 cycles/deg
    'orientationDegs', 45, ...              % 45 degrees
    'widthDegs', 0.4, ...                   % 0.x4 x 0.4 size
    'contrast', 0.8,...                     % 80% Michelson contrast
    'meanLuminanceCdPerM2', 40);            % 40 cd/m2 mean luminance

%% Generate an ISETBio scene describing this stimulus
scene = generateStimulusScene(stimParams, presentationDisplay);

%% Realize this scene into the particular LCD display
realizedScene = realizeSceneInDisplay(scene, presentationDisplay);

%% Generate wavefront-aberration derived human optics
opticalImage = oiCreate('wvf human');

%% Compute the retinal image
opticalImage = oiCompute(opticalImage, realizedScene);

%% Generate a hexagonal cone mosaic with ecc-based cone quantal efficiency
coneMosaic = coneMosaicHex(7, ...             % hex lattice sampling factor
   'fovDegs', stimParams.widthDegs, ...       % match mosaic width to stimulus size
   'eccBasedConeDensity', true, ...           % cone density varies with eccentricity
   'eccBasedConeQuantalEfficiency', true, ... % cone quantal efficiency varies with eccentricity
   'integrationTime', 10/1000, ...            % 10 msec integration time
   'maxGridAdjustmentIterations', 50);        % terminate iterative lattice adjustment after 50 iterations

%% Compute cone excitations in response to the stimulus 
nInstances = 3;   % generate 3 response instances
coneExcitations = coneMosaic.compute(opticalImage, 'emPath', zeros(nInstances, 1, 2));

%% Step 9. Display components
% Decide which wavelengths to visualize
displayedWavelengths = [450:25:725];         % wavelength bands to visualize
% Get range of scene photons
scenePhotons = sceneGet(scene, 'photons');
realizedScenePhotons = sceneGet(realizedScene, 'photons');
maxPhotons = max([ max(scenePhotons(:)) max(realizedScenePhotons(:))]);

%% Visualize scene 
figNo = 1;
visualizeScene(scene, maxPhotons, displayedWavelengths , figNo, 'scene');

%% Visualize display
figNo = figNo + 1;
visualizeDisplay(presentationDisplay, figNo);

%% Visualize scene rendered on display
figNo = figNo + 1;
visualizeScene(realizedScene, maxPhotons, displayedWavelengths, figNo, 'realized scene');

%% Compare the 2 scenes
figNo = figNo + 1;
visualizeSpectralSlices(scene, realizedScene, maxPhotons, displayedWavelengths, figNo);

%% Visualize the PSF
figNo = figNo + 1;
visualizePSF(opticalImage, displayedWavelengths, figNo);

%% Visualize the optical image
figNo = figNo + 1;
visualizeOpticalImage(opticalImage, displayedWavelengths, figNo);

%% Visualize cone mosaic and cone excitation response
figNo = figNo + 1;
visualizeResponses(coneMosaic, coneExcitations, figNo);
 