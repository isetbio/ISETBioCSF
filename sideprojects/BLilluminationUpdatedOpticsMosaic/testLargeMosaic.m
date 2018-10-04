function testLargeMosaic

    ibioDataDir = '/Volumes/DropBoxDisk/Dropbox/Dropbox (Aguirre-Brainard Lab)/IBIO_data';
    %ibioDataDir = '/Volumes/SamsungT3/Dropbox/AguirreBrainardLabsDropbox/IBIO_data';
    %ibioDataDir = '/Volumes/IthakasPassport/DropboxLab/AquirreBrainardLab/IBIO_data'

    sceneRootDir = fullfile(ibioDataDir,'BLIlluminationDiscrimination/SceneData');
    oiRootDir = fullfile(ibioDataDir, 'BLIlluminationDiscrimination/OpticalImageDataNoBlur');
    
    %opticsModel = 'ThibosDefaultSubject3MMPupil';
    opticsModel = 'DeltaFunction';
    
    sceneNames = {...
        'Constant_CorrectSize'...
      %  'NM1_CorrectSize' ...		
      %  'NM2_CorrectSize' ...		
      %  'Neutral_CorrectSize' ...
        };
    illuminationNames = { ...
        'BlueIllumination' ...
        'GreenIllumination' ...	
        'RedIllumination' ...		
        'Standard' ...		
        'YellowIllumination' ...
        };
    %illuminationNames = {illuminationNames{4}};
    
    sceneNamesNum = 1; % numel(sceneNames)
    illuminationsNum = 1; % numel(illuminationNames)
    for sceneIndex = 1:sceneNamesNum
        sceneName = sceneNames{sceneIndex};
        for illumIndex = 1:illuminationsNum
        
            illuminationName = illuminationNames{illumIndex};
        
            sceneDir = fullfile(sceneRootDir, sceneName, illuminationName);
            oiDir = fullfile(oiRootDir, sceneName, illuminationName);
        
            listings = dir(sprintf('%s/*.mat', sceneDir));
        
            readOneOI = false;
            
            for k = 1:numel(listings)
                if (contains(listings(k).name, '.mat'))  && (~readOneOI)
                    % Source filename
                    
                    oiFileName = fullfile(oiDir, listings(k).name);
                    load(oiFileName, 'oi');
                    fprintf('oisize: %f x %f degs\n', oiGet(oi,'hfov'), oiGet(oi, 'vfov'));
                    readOneOI = true;
                end
            end % k
        end % illumIndex
    end % sceneIndex
    
    resamplingFactor = 5;
    integrationTime = 50/1000;
    
    if (1==2)
    cReg = coneMosaicHex(resamplingFactor, ...
        'fovDegs', [oiGet(oi,'hfov') oiGet(oi, 'vfov')]/10, ...
        'integrationTime', integrationTime);
    fprintf('Cones num: %d\n', numel(find(cReg.pattern > 1)));
    
    cReg.compute(oi);
    pause
    end
    
    cEccBased = coneMosaicHex(resamplingFactor, ...
        'fovDegs', [oiGet(oi,'hfov') oiGet(oi, 'vfov')], ...
        'integrationTime', integrationTime, ...
        'eccBasedConeDensity', true, ...
        'eccBasedConeEfficiency', true, ...
        'maxGridAdjustmentIterations', 3);
    fprintf('Cones num: %d\n', numel(find(cEccBased.pattern > 1)));
    
    save('largeMosaic.mat', 'cEccBased');
end
