function testConeContrasts = testConeContrastsFromTestDirectionParams(rParams,testDirectionParams)
% testConeContrasts = testConeContrastsFromTestDirectionParams(rParams,testDirectionParams)
%
% Convert various test direction specifications to a set of test cone
% contrasts.  
%
% Scales contrasts to lie within monitor gamut, which provides a better
% scale for guessing contrast multipliers than just normalizing to unit
% length.

%% Define test cone contrasts
%
% Chromatic directions in L/M plane.  It's a little easier to think in
% terms of angles.  
%
% Then find highest within gamut vector for each direction
switch (testDirectionParams.type)
    case 'LMPlaneInstance'
        LMangles = (0:testDirectionParams.deltaAngle:180-testDirectionParams.deltaAngle)/180*pi;
        for angleIndex = 1:numel(LMangles)
            theta = LMangles(angleIndex);
            baseTestConeContrastDirs(:,angleIndex) = testDirectionParams.baseStimulusLength*[cos(theta) sin(theta) 0.0]';
            
            % Find the highest in gamut cone contrast and define cone contrast
            % vector to be just under this length.
            colorModulationParamsTemp = rParams.colorModulationParams;
            colorModulationParamsTemp.coneContrasts = baseTestConeContrastDirs(:,angleIndex);
            colorModulationParamsTemp.contrast = 1;
            [~,contrastScaleFactor(angleIndex)] = colorGaborSceneCreate(rParams.gaborParams,colorModulationParamsTemp,true);
            testConeContrasts(:,angleIndex) = 0.98*contrastScaleFactor(angleIndex)*baseTestConeContrastDirs(:,angleIndex);
        end
end
