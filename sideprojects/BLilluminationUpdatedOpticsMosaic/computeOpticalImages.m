function computeOpticalImages

%     visualizePSF = true;
%     horizontalFOV = 16;
%     pupilDiamMM = 6;
%     generateOI(horizontalFOV, pupilDiamMM, visualizePSF);
%     return;
    
    sceneRootDir = '/Volumes/DropBoxDisk/Dropbox/Dropbox (Aguirre-Brainard Lab)/IBIO_data/BLIlluminationDiscrimination/SceneData';
    oiRootDir = '/Volumes/DropBoxDisk/Dropbox/Dropbox (Aguirre-Brainard Lab)/IBIO_data/BLIlluminationDiscrimination/OpticalImageData';
    
    sceneNames = {...
        'Constant_CorrectSize_Blue'...
        %'NM1_CorrectSize' ...		
        %'NM2_CorrectSize' ...		
        %'Neutral_CorrectSize' ...
        };
    illuminationNames = { ...
        'BlueIllumination' ...
       % 'GreenIllumination' ...	
       % 'RedIllumination' ...		
       % 'Standard' ...		
       % 'YellowIllumination' ...
        };
    
    localDir = pwd;
    for sceneIndex = 1:numel(sceneNames)
        sceneName = sceneNames{sceneIndex};
        if (~isdir(fullfile(oiRootDir, sceneName)))
            cd(oiRootDir)
            mkdir(sceneName);
            cd(localDir);
        end
        if (strcmp(sceneName, 'Constant_CorrectSize'))
            % 6 mm for modeled experiment 1
            pupilDiamMM = 6;
        else
            % 5 mm for modeled experiment 2
            pupilDiamMM = 5;
        end
        for illumIndex = 1:numel(illuminationNames)
            illuminationName = illuminationNames{illumIndex};
            if (~isdir(fullfile(oiRootDir, sceneName, illuminationName)))
                cd(fullfile(oiRootDir, sceneName));
                mkdir(illuminationName);
                cd(localDir);
            end
        
            fprintf('Using %d mm pupil for scenes named: %s/%s\n', pupilDiamMM,sceneName, illuminationName);
            sceneDir = fullfile(sceneRootDir, sceneName, illuminationName);
            oiDir = fullfile(oiRootDir, sceneName, illuminationName);
            
            listings = dir(sprintf('%s/*.mat', sceneDir));
            for k = 1:numel(listings)
                if (contains(listings(k).name, '.mat'))  
                    % Source filename
                    sceneFileName = fullfile(sceneDir, listings(k).name);
                    % Destination filename
                    oiFileName  = fullfile(oiDir, listings(k).name);
                    % Load the scene
                    load(sceneFileName, 'scene');
                    theScene = scene;
                    clear 'scene';
                    horizontalFOV = sceneGet(theScene, 'hfov');
                    visualizePSF = false;
                    theOI = generateOI(horizontalFOV, pupilDiamMM, visualizePSF);
                    oi = oiCompute(theOI, theScene);
                    save(oiFileName, 'oi');
                    fprintf('Saved oi in %s\n', oiFileName);
                    visualizeSceneAndOI = false;
                    if (visualizeSceneAndOI)
                        figure(10);
                        subplot(1,2,1)
                        imagesc(sceneGet(theScene, 'rgb image'));
                        axis 'image'
                        subplot(1,2,2)
                        imagesc(oiGet(oi, 'rgb image'));
                        axis 'image'
                    end
                end
            end
        end % illumIndex
    end % for sceneIndex
end

function theOI = generateOI(horizontalFOV, pupilDiamMM, visualizePSF)

    % Generate human optics
    oiParams = struct(...
        'opticsModel', 'ThibosDefaultSubject3MMPupil', ...
        'wavefrontSpatialSamples', 261*2+1, ...
        'pupilDiamMm', pupilDiamMM, ...
        'umPerDegree', 300);
    theOI = oiWithCustomOptics(oiParams.opticsModel, oiParams.wavefrontSpatialSamples, oiParams.pupilDiamMm, oiParams.umPerDegree);

    % Set the FOV
    theOI = oiSet(theOI,'h fov',horizontalFOV);

    % Set the fNumber
    focalLength = oiGet(theOI,'distance');
    desiredFNumber = focalLength/(oiParams.pupilDiamMm/1000);
    theOI  = oiSet(theOI ,'optics fnumber',desiredFNumber);
    
    if (visualizePSF)
    visualizePSFfromOI(theOI, oiParams.umPerDegree, ...
                'colormapToUse', gray(1024), ...
                'visualizedWavelengths', [430 490 550 610 670], ...
                'rows', 1, 'cols', 5, ...
                'labelLastPSF', true, ...
                'displayWavelengthInTitle', true);
    end   
end


