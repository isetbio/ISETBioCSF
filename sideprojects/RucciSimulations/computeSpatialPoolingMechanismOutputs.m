function computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, parforWorkers, resourcesDir)
    
    % Load responses to null (zero contrast) stimulus
    fname = fullfile(resourcesDir, sprintf('%s_nTrials_%d.mat', 'zeroContrast',  nTrials));
    load(fname, 'coneExcitations', 'photoCurrents', 'emPathsDegs', 'timeAxis');

    % Get the last time point at the mean response
    nullStimulusMeanConeExcitations = coneExcitations(1,:,end);
    nullStimulusMeanPhotocurrents = photoCurrents(1,:,end);
    
    if (contains(stimDescriptor, 'highFrequency'))
        poolingKernelLinear     = spatialPoolingKernels.highFrequencyPoolingWeightsLinear;
        poolingKernelQuadrature = spatialPoolingKernels.highFrequencyPoolingWeightsQuadrature;
        poolingKernelOrthoLinear     = spatialPoolingKernels.highFrequencyOrthoPoolingWeightsLinear;
        poolingKernelOrthoQuadrature = spatialPoolingKernels.highFrequencyOrthoPoolingWeightsQuadrature;
        
    elseif (contains(stimDescriptor, 'lowFrequency'))
        poolingKernelLinear     = spatialPoolingKernels.lowFrequencyPoolingWeightsLinear;
        poolingKernelQuadrature = spatialPoolingKernels.lowFrequencyPoolingWeightsQuadrature;
        poolingKernelOrthoLinear     = spatialPoolingKernels.lowFrequencyOrthoPoolingWeightsLinear;
        poolingKernelOrthoQuadrature = spatialPoolingKernels.lowFrequencyOrthoPoolingWeightsQuadrature;
    else
        error('Unknown stimDescriptor: ''%s''.', stimDescriptor);
    end
    
    % Unit volume
    poolingKernelLinear = poolingKernelLinear / sum(abs(poolingKernelLinear));
    poolingKernelQuadrature = poolingKernelQuadrature / sum(abs(poolingKernelQuadrature));
    poolingKernelOrthoLinear = poolingKernelOrthoLinear / sum(abs(poolingKernelOrthoLinear));
    poolingKernelOrthoQuadrature = poolingKernelOrthoQuadrature / sum(abs(poolingKernelOrthoQuadrature));
    
    % Compute modulated mosaic responses for all other OIs and all contrast levels
    nContrasts = numel(contrastLevels);
    
    % Preallocate memory for results
    energyConeExcitationResponse.output = zeros(nContrasts, nTrials, size(coneExcitations,3));
    energyConeExcitationResponse.orthoOutput = zeros(nContrasts, nTrials, size(coneExcitations,3));
    energyPhotoCurrentResponse.output = zeros(nContrasts, nTrials, size(photoCurrents,3));
    energyPhotoCurrentResponse.orthoOutput = zeros(nContrasts, nTrials, size(photoCurrents,3));
    
    for theContrastLevel = 1:nContrasts
        % Load mosaic responses
        fname = fullfile(resourcesDir, sprintf('%s_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', stimDescriptor, contrastLevels(theContrastLevel), analyzedNoiseInstance, nTrials));
        load(fname, 'coneExcitations', 'photoCurrents');
        
        % Compute cone excitation energy response for the standard orientation filter
        energyConeExcitationResponse.output(theContrastLevel,:,:) = ...
            computeEnergyResponse(coneExcitations, nullStimulusMeanConeExcitations, ...
            poolingKernelLinear, poolingKernelQuadrature);
        
        % Compute cone excitation energy response for the orthogonal orientation filter
        energyConeExcitationResponse.orthoOutput(theContrastLevel,:,:) = ...
            computeEnergyResponse(coneExcitations, nullStimulusMeanConeExcitations, ...
            poolingKernelOrthoLinear, poolingKernelOrthoQuadrature);
        
        % Compute photocurrent energy response for the standard orientation filter
        energyPhotoCurrentResponse.output(theContrastLevel,:,:) = ...
            computeEnergyResponse(photoCurrents, nullStimulusMeanPhotocurrents, ...
            poolingKernelLinear, poolingKernelQuadrature);
        
        % Compute photocurrent energy response for the orthogonal orientation filter
        energyPhotoCurrentResponse.orthoOutput(theContrastLevel,:,:) = ...
            computeEnergyResponse(photoCurrents, nullStimulusMeanPhotocurrents, ...
            poolingKernelOrthoLinear, poolingKernelOrthoQuadrature);  
    end
    
    fname = fullfile(resourcesDir, sprintf('energyResponse_%s_instance_%1.0f_nTrials_%d.mat', stimDescriptor, analyzedNoiseInstance, nTrials));
    fprintf('Saving energy responses from %d trials for %s stimulus to %s\n', nTrials, stimDescriptor, fname);
    save(fname, 'energyConeExcitationResponse', 'energyPhotoCurrentResponse', 'emPathsDegs', 'timeAxis', '-v7.3');

end


function energyResponse = computeEnergyResponse(response, nullStimulusResponse, poolingKernelLinear, poolingKernelQuadrature)
    % Compute modulated mosaic responses by subtracting the MEAN response to the null stimulus
    response = bsxfun(@minus, response, nullStimulusResponse);

    % Comptute dot product along all cones for the quadrature pair of pooling kernels
    subunit1Response = squeeze(sum(bsxfun(@times, poolingKernelLinear, response),2));
    subunit2Response = squeeze(sum(bsxfun(@times, poolingKernelQuadrature, response),2));

    % Sum squared outputs
    energyResponse = subunit1Response.^2 + subunit2Response.^2;
end