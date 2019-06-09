function fName = oisDataFileName(meanLuminanceCdPerM2, contrastLevels,resourcesDir)
    minC = min(contrastLevels);
    maxC = max(contrastLevels);
    nContrasts = numel(contrastLevels);
    fName = fullfile(resourcesDir, sprintf('ois_luminance_%2.1f_minC_%2.4f_maxC_%2.4f_nC_%d.mat', meanLuminanceCdPerM2, minC, maxC, nContrasts));
end

