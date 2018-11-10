%% Set the viewing distance
viewingDistanceMeters = 0.57;

%% Step 1. Generate a display for presenting stimuli and place it at the desired viewing distance
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', viewingDistanceMeters);

%% Step 2. Specify stimulus params
stimParams = struct(...
    'spatialFrequencyCyclesPerDeg', 10, ... 
    'orientationDegs', 45, ...
    'widthDegs', 0.4, ...
    'contrast', 0.8, ...
    'meanLuminanceCdPerM2', 40);

%% Step 3. Generate a scene describing the stimulus
scene = generateStimulusScene(stimParams, presentationDisplay);

%% Step 4. Realize this scene into a particular display
realizedScene = realizeSceneInDisplay(scene, presentationDisplay);

%% Step 5. Generate human optics
opticalImage = oiCreate('wvf human');

%% Step 6. Compute retinal image
opticalImage = oiCompute(opticalImage, realizedScene);

%% Step 7. Generate hexagonal cone mosaic
coneMosaic = coneMosaicHex(7, ...
    'fovDegs', stimParams.widthDegs, ...
    'eccBasedConeDensity', true, ...
    'eccBasedConeQuantalEfficiency', true, ...
    'maxGridAdjustmentIterations', 50);

%% Compute using a 10 ms integration time
coneMosaic.integrationTime = 10/1000;

%% Step 8. Compute cone excitations
nInstances = 3;
coneExcitations = coneMosaic.compute(opticalImage, 'emPath', zeros(nInstances, 1, 2));

%% Step 9. Display components
displayedWavelengths = [450:25:725];
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
%figNo = figNo + 1;
%visualizeSpectralSlices(scene, realizedScene, maxPhotons, displayedWavelengths, figNo);

%% Visualize the PSF
figNo = figNo + 1;
visualizePSF(opticalImage, displayedWavelengths, figNo);

%% Visualize the optical image
figNo = figNo + 1;
visualizeOpticalImage(opticalImage, displayedWavelengths, figNo);

%% Visualize cone mosaic and cone excitation response
figNo = figNo + 1;
visualizeResponses(coneMosaic, coneExcitations, figNo);
 