% colorGaborSceneCreateTest
%
% Test the routine that creats a color gabor scene.
%
% 7/7/16  dhb  Wrote it.

%% Clear
clear; close all;

%% Add project toolbox to Matlab path
AddToMatlabPathDynamically(fullfile(fileparts(which(mfilename)),'../toolbox')); 

%% Make sure we can make a basic scene.
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
gaborParams.monitorFile = 'CRT-HP';
gaborParams.viewingDistance = 1;
gaborScene = colorGaborSceneCreate(gaborParams);
vcAddAndSelectObject(gaborScene);sceneWindow;