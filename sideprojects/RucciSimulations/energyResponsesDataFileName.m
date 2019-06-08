function fName = energyResponsesDataFileName(stimDescriptor,  analyzedNoiseInstance, nTrials, eyePosition, resourcesDir)
    fName = fullfile(resourcesDir, sprintf('energyResponse_%s_instance_%1.0f_nTrials_%d.mat', stimDescriptor, analyzedNoiseInstance, nTrials));
    if (strcmp(eyePosition, 'stabilized'))
        fName = fullfile(resourcesDir, sprintf('energyResponse_%s_instance_%1.0f_nTrials_%d_STABILIZED.mat', stimDescriptor, analyzedNoiseInstance, nTrials));
    end
end

