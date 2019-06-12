function fName = performanceDataFileName(stimDescriptor, analyzedNoiseInstance, nTrials, eyePosition, resourcesDir)
    fName = fullfile(resourcesDir, sprintf('performance_%s_instance_%1.0f_nTrials_%d_%s.mat', stimDescriptor, analyzedNoiseInstance, nTrials, upper(eyePosition)));
end
