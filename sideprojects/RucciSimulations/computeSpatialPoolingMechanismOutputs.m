function computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, parforWorkers, resourcesDir)
    
    % Load responses to null (zero contrast) stimulus
    fname = fullfile(resourcesDir, sprintf('%s_nTrials_%d.mat', 'zeroContrast',  nTrials));
    load(fname, 'coneExcitations', 'emPathsDegs', 'timeAxis');

    % Get the last time point at the mean response
    nullStimulusMeanConeExcitations = coneExcitations(1,:,end);
    
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
    %poolingKernelLinear = poolingKernelLinear / sum(abs(poolingKernelLinear));
    %poolingKernelQuadrature = poolingKernelQuadrature / sum(abs(poolingKernelQuadrature));
    %poolingKernelOrthoLinear = poolingKernelOrthoLinear / sum(abs(poolingKernelOrthoLinear));
    %poolingKernelOrthoQuadrature = poolingKernelOrthoQuadrature / sum(abs(poolingKernelOrthoQuadrature));
    
    % Compute modulated mosaic responses for all other OIs and all contrast levels
    nContrasts = numel(contrastLevels);
    
    % Preallocate memory for results
    energyResponse.output = zeros(nContrasts, nTrials, size(coneExcitations,3));
    energyResponse.orthoOutput = zeros(nContrasts, nTrials, size(coneExcitations,3));
    
    figure(9876); clf;
    spatialSupportDegs = 0.5;
    ax = subplot(2,2,1);
    visualizePoolingKernel(ax, spatialPoolingKernels.coneLocsDegs, ...
        spatialPoolingKernels.coneAperture, poolingKernelLinear, spatialSupportDegs, true)
    
    ax = subplot(2,2,2);
    visualizePoolingKernel(ax, spatialPoolingKernels.coneLocsDegs, ...
        spatialPoolingKernels.coneAperture, poolingKernelOrthoLinear, spatialSupportDegs, true)
    
    stimDescriptor
    pause
    
    for theContrastLevel = 1:1 % nContrasts
        % Load mosaic responses
        fname = fullfile(resourcesDir, sprintf('%s_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', stimDescriptor, contrastLevels(theContrastLevel), analyzedNoiseInstance, nTrials));
        load(fname, 'coneExcitations');
        
        % Compute modulated mosaic responses by subtracting the MEAN response to the null stimulus
        modulatedConeExcitations = bsxfun(@minus, coneExcitations, nullStimulusMeanConeExcitations);
        

        
        for timeBin = 1:80
        ax = subplot(2,2,3)
        visualizePoolingKernel(ax, spatialPoolingKernels.coneLocsDegs, ...
            spatialPoolingKernels.coneAperture, squeeze(modulatedConeExcitations(1,:,timeBin)), spatialSupportDegs, true)
        title(sprintf('time: %d of %d', timeBin,80));
        
        ax = subplot(2,2,4)
        visualizePoolingKernel(ax, spatialPoolingKernels.coneLocsDegs, ...
            spatialPoolingKernels.coneAperture, squeeze(coneExcitations(1,:,timeBin)), spatialSupportDegs, true)
        drawnow;
        pause(0.1);
        end
    
    
        % Comptute dot product along all cones (standard orientation filter)
        subunit1Responses = squeeze(sum(bsxfun(@times, poolingKernelLinear, modulatedConeExcitations),2));
        subunit2Responses = squeeze(sum(bsxfun(@times, poolingKernelQuadrature, modulatedConeExcitations),2));
        
        % Compute energy (standard orientation filter)
        energyResponse.output(theContrastLevel,:,:) = subunit1Responses.^2 + subunit2Responses.^2;
        
        % Comptute dot product along all cones (orthogonal orientation filter)
        subunit1Responses = squeeze(sum(bsxfun(@times, poolingKernelOrthoLinear, modulatedConeExcitations),2));
        subunit2Responses = squeeze(sum(bsxfun(@times, poolingKernelOrthoQuadrature, modulatedConeExcitations),2));
        
        % Compute energy (orthogonal orientation filter)
        energyResponse.orthoOutput(theContrastLevel,:,:) = subunit1Responses.^2 + subunit2Responses.^2;
    end
    
    fname = fullfile(resourcesDir, sprintf('energyResponse_%s_instance_%1.0f_nTrials_%d.mat', stimDescriptor, analyzedNoiseInstance, nTrials));
    fprintf('Saving energy responses from %d trials for %s stimulus to %s\n', nTrials, stimDescriptor, fname);
    save(fname, 'energyResponse', 'emPathsDegs', 'timeAxis', '-v7.3');

end


