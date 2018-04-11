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
switch (testDirectionParams.instanceType)
    case 'LMPlane'
        baseTestConeContrastDirs = zeros(3, testDirectionParams.nAngles);
        for angleIndex = 1:testDirectionParams.nAngles
            theta = (pi/180)*(testDirectionParams.startAngle + (angleIndex-1)*testDirectionParams.deltaAngle);
            baseTestConeContrastDirs(:,angleIndex) = testDirectionParams.baseStimulusLength*[cos(theta) sin(theta) 0.0]';
        end
        
    case 'LMSPlane' 
        % Generate meshgrid from all tested azimuth and elevation angles
        azimuthAngles = testDirectionParams.startAzimuthAngle + (0:(testDirectionParams.nAzimuthAngles-1))*testDirectionParams.deltaAzimuthAngle;
        elevationAngles = testDirectionParams.startElevationAngle + (0:(testDirectionParams.nElevationAngles-1))*testDirectionParams.deltaElevationAngle;
        [azimuths,elevations] = meshgrid(azimuthAngles,elevationAngles);
   
        baseTestConeContrastDirs = zeros(3, numel(azimuths));
        for angleIndex = 1:numel(azimuths)
            azAngle = azimuths(angleIndex);
            elAngle = elevations(angleIndex);
            % Compute cone contrasts from azimuth and elevation
            [cL, cM, cS] = sph2cart(azAngle/180*pi, elAngle/180*pi, 1);
            baseTestConeContrastDirs(:,angleIndex) = testDirectionParams.baseStimulusLength * [cL, cM, cS]';
        end
end

for angleIndex = 1:size(baseTestConeContrastDirs,2)
    % Find the highest in gamut cone contrast and define cone contrast
    % vector to be just under this length.
    colorModulationParamsTemp = rParams.colorModulationParams;
    colorModulationParamsTemp.coneContrasts = baseTestConeContrastDirs(:,angleIndex);
    colorModulationParamsTemp.contrast = 1;
    [~,contrastScaleFactor(angleIndex)] = colorSceneCreate(rParams.spatialParams,rParams.backgroundParams,colorModulationParamsTemp,rParams.oiParams,true);
    testConeContrasts(:,angleIndex) = 0.98*contrastScaleFactor(angleIndex)*baseTestConeContrastDirs(:,angleIndex);
    fprintf('Max in-gamut cone contrast for direction #%d: %f %f %f\n', angleIndex, testConeContrasts(1,angleIndex), testConeContrasts(2,angleIndex), testConeContrasts(3,angleIndex));
end
end

 
            
