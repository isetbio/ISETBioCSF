%% Visualize results at selected wavelengths
displayedWavelengths = [450:25:725];
scenePhotons = sceneGet(scene, 'photons');
realizedScenePhotons = sceneGet(realizedScene, 'photons');
maxPhotons = max([ max(scenePhotons(:)) max(realizedScenePhotons(:))]);

%% Visualize components of the desired scene
figNo = 1;
visualizeScene(scene, maxPhotons, displayedWavelengths , figNo, 'scene');

%% Visualize components of the LCD display
figNo = figNo + 1;
visualizeDisplay(presentationDisplay, figNo);

%% Visualize components of the LCD display - rendered scene
figNo = figNo + 1;
visualizeScene(realizedScene, maxPhotons, displayedWavelengths, figNo, 'realized scene');

%% Visualize the PSF
figNo = figNo + 1;
visualizePSF(opticalImage, displayedWavelengths, figNo);

%% Visualize the optical image
figNo = figNo + 1;
visualizeOpticalImage(opticalImage, displayedWavelengths, figNo);

%% Visualize the cone mosaic and the cone excitation responses
figNo = figNo + 1;
visualizeResponses(coneMosaic, coneExcitations, figNo);
 