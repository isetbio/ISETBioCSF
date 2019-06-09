function computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, eyePosition, parforWorkers, resourcesDir)
    
    % Load responses to null (zero contrast) stimulus
    fName = fullfile(resourcesDir, sprintf('%s_nTrials_%d.mat', 'zeroContrast',  nTrials));
    load(fName, 'coneExcitations', 'photoCurrents', 'emPathsDegs', 'timeAxis');

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
    energyConeExcitationResponseOutput = zeros(nContrasts, nTrials, size(coneExcitations,3));
    energyConeExcitationResponseOrthoOutput = zeros(nContrasts, nTrials, size(coneExcitations,3));
    energyPhotoCurrentResponseOutput = zeros(nContrasts, nTrials, size(photoCurrents,3));
    energyPhotoCurrentResponseOrthoOutput = zeros(nContrasts, nTrials, size(photoCurrents,3));
    
    parfor theContrastLevel = 1:nContrasts
        fprintf('Computing energy response for contrast %d\n', theContrastLevel);
        % Load mosaic responses
        [coneExcitations, photoCurrents] = ...
            retrieveData(stimDescriptor, contrastLevels(theContrastLevel), analyzedNoiseInstance, nTrials, eyePosition, resourcesDir);
        
        % Compute cone excitation energy response for the standard orientation filter
        energyConeExcitationResponseOutput(theContrastLevel,:,:) = ...
            computeEnergyResponse(coneExcitations, nullStimulusMeanConeExcitations, ...
            poolingKernelLinear, poolingKernelQuadrature);
        
        % Compute cone excitation energy response for the orthogonal orientation filter
        energyConeExcitationResponseOrthoOutput(theContrastLevel,:,:) = ...
            computeEnergyResponse(coneExcitations, nullStimulusMeanConeExcitations, ...
            poolingKernelOrthoLinear, poolingKernelOrthoQuadrature);
        
        % Compute photocurrent energy response for the standard orientation filter
        energyPhotoCurrentResponseOutput(theContrastLevel,:,:) = ...
            computeEnergyResponse(photoCurrents, nullStimulusMeanPhotocurrents, ...
            poolingKernelLinear, poolingKernelQuadrature);
        
        % Compute photocurrent energy response for the orthogonal orientation filter
        energyPhotoCurrentResponseOrthoOutput(theContrastLevel,:,:) = ...
            computeEnergyResponse(photoCurrents, nullStimulusMeanPhotocurrents, ...
            poolingKernelOrthoLinear, poolingKernelOrthoQuadrature);  
    end
    
    fName = energyResponsesDataFileName(stimDescriptor, analyzedNoiseInstance, nTrials, eyePosition, resourcesDir);
    fprintf('Saving energy responses from %d trials for %s stimulus to %s\n', nTrials, stimDescriptor, fName);
    save(fName, 'energyConeExcitationResponseOutput', 'energyConeExcitationResponseOrthoOutput', ...
                'energyPhotoCurrentResponseOutput', 'energyPhotoCurrentResponseOrthoOutput', ...
                'emPathsDegs', 'timeAxis', '-v7.3');

end

function [coneExcitations, photoCurrents] = retrieveData(stimDescriptor, theContrastLevel, analyzedNoiseInstance, nTrials, eyePosition, resourcesDir)
    fName = coneMosaicResponsesDataFileName(stimDescriptor, theContrastLevel, analyzedNoiseInstance, nTrials, eyePosition, resourcesDir);
    load(fName, 'coneExcitations', 'photoCurrents'); 
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