function fName = coneMosaicResponsesDataFileName(stimDescriptor, contrastLevel, analyzedNoiseInstance, nTrials, eyePosition, resourcesDir)
    fName = fullfile(resourcesDir, sprintf('%s_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', stimDescriptor, contrastLevel, analyzedNoiseInstance, nTrials));
    if (strcmp(eyePosition, 'stabilized'))
        fName = fullfile(resourcesDir, sprintf('%s_contrast_%2.4f_instance_%1.0f_nTrials_%d_STABILIZED.mat', stimDescriptor, contrastLevel, analyzedNoiseInstance, nTrials));
    end
end

