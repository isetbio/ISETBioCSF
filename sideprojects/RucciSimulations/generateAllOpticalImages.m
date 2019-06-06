    function generateAllOpticalImages(nullScene, lowFrequencyScenes, highFrequencyScenes, lowFrequencyScenesOrtho, highFrequencyScenesOrtho, contrastLevels, noiseInstances, meanLuminanceCdPerM2, resourcesDir)
    nContrasts = numel(contrastLevels);
    theOI = oiCreate('wvf human');
    
    % The zero contrast scene oi
    nullSceneOI = oiCompute(theOI, nullScene);
        
    % The other ois
    for theContrastLevel = 1:nContrasts
        for theInstance = 1:noiseInstances       
            lowFrequencyOIs{theContrastLevel, theInstance} = oiCompute(theOI, ...
            lowFrequencyScenes{theContrastLevel, theInstance});
        
            lowFrequencyOIsOrtho{theContrastLevel, theInstance} = oiCompute(theOI, ...
            lowFrequencyScenesOrtho{theContrastLevel, theInstance});
        
            highFrequencyOIs{theContrastLevel, theInstance} = oiCompute(theOI, ...
            highFrequencyScenes{theContrastLevel, theInstance});
        
            highFrequencyOIsOrtho{theContrastLevel, theInstance} = oiCompute(theOI, ...
            highFrequencyScenesOrtho{theContrastLevel, theInstance});
        end
    end
    
    fName = fullfile(resourcesDir, sprintf('ois_luminance_%2.1f.mat', meanLuminanceCdPerM2));
    save(fName, 'nullSceneOI', 'lowFrequencyOIs', 'highFrequencyOIs', ...
         'lowFrequencyOIsOrtho', 'highFrequencyOIsOrtho', 'contrastLevels', 'noiseInstances', '-v7.3');        
end

