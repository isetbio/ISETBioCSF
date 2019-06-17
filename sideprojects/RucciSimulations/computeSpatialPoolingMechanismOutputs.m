function computeSpatialPoolingMechanismOutputs(spatialPoolingKernels, stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, eyePosition, parforWorkers, resourcesDir)
    
    % Load responses to null (zero contrast) stimulus
    fName = fullfile(resourcesDir, sprintf('%s_nTrials_%d.mat', 'zeroContrast',  nTrials));
    load(fName, 'coneExcitations', 'photoCurrents', 'emPathsDegs', 'timeAxis');

    temporalKernels = computeTemporalKernel(timeAxis);
    
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
    
    for theContrastLevel = 1:nContrasts
        fprintf('Computing energy response for contrast %d\n', theContrastLevel);
        % Load mosaic responses
        [coneExcitations, photoCurrents] = ...
            retrieveData(stimDescriptor, contrastLevels(theContrastLevel), analyzedNoiseInstance, nTrials, eyePosition, resourcesDir);
        
        % Compute cone excitation energy response for the standard orientation filter
        energyConeExcitationResponseOutput(theContrastLevel,:,:) = ...
            computeEnergyResponse(coneExcitations, nullStimulusMeanConeExcitations, ...
            poolingKernelLinear, poolingKernelQuadrature, temporalKernels);
        
        % Compute cone excitation energy response for the orthogonal orientation filter
        energyConeExcitationResponseOrthoOutput(theContrastLevel,:,:) = ...
            computeEnergyResponse(coneExcitations, nullStimulusMeanConeExcitations, ...
            poolingKernelOrthoLinear, poolingKernelOrthoQuadrature, temporalKernels);
        
        % Compute photocurrent energy response for the standard orientation filter
        energyPhotoCurrentResponseOutput(theContrastLevel,:,:) = ...
            computeEnergyResponse(photoCurrents, nullStimulusMeanPhotocurrents, ...
            poolingKernelLinear, poolingKernelQuadrature, temporalKernels);
        
        % Compute photocurrent energy response for the orthogonal orientation filter
        energyPhotoCurrentResponseOrthoOutput(theContrastLevel,:,:) = ...
            computeEnergyResponse(photoCurrents, nullStimulusMeanPhotocurrents, ...
            poolingKernelOrthoLinear, poolingKernelOrthoQuadrature, temporalKernels);  
    end
    
    fName = energyResponsesDataFileName(stimDescriptor, analyzedNoiseInstance, nTrials, eyePosition, resourcesDir);
    fprintf('Saving energy responses from %d trials for %s stimulus to %s\n', nTrials, stimDescriptor, fName);
    save(fName, 'energyConeExcitationResponseOutput', 'energyConeExcitationResponseOrthoOutput', ...
                'energyPhotoCurrentResponseOutput', 'energyPhotoCurrentResponseOrthoOutput', ...
                'emPathsDegs', 'timeAxis', '-v7.3');

end

function temporalKernels = computeTemporalKernel(timeAxis)
    k = 50;
    dt = timeAxis(2)-timeAxis(1);
    temporalSupport = 0:dt:(400/1000);
    kT = k * temporalSupport;
    N = 3;
    temporalKernels(1,:) = (kT.^N) .* (exp(-kT)) .* (1/factorial(N) - (kT.^2)/factorial(N+2));
    N = 5;
    temporalKernels(2,:) = (kT.^N) .* (exp(-kT)) .* (1/factorial(N) - (kT.^2)/factorial(N+2));
    
    ft = abs(fftshift(fft(squeeze(temporalKernels(1,:)),4096)));
    maxTF = 1 / (2*dt);
    temporalFreqSupport = ((1:numel(ft))/numel(ft) -0.5)*2 * maxTF;
        
    xticks = [0:100:500];
    xTickLabels = sprintf('%2.0f\n', xticks);
    figure()
    subplot(1,2,1)
    plot(temporalSupport*1000, temporalKernels(1,:), 'r-', 'LineWidth', 1.5); hold on
    plot(temporalSupport*1009, temporalKernels(2,:), 'b-', 'LineWidth', 1.5);
    
    set(gca, 'YTickLabel', {}, 'XTick', xticks, 'XTickLabel', xTickLabels, 'XLim', [0 400], 'FontSize', 14);
    xlabel('\it time (msec)');
    ylabel('');
    axis 'square';
    grid on; box on
    
    xticks = [0.3 1 3 10 30 100];
    xTickLabels = cell(1, numel(xticks));
    for k = 1:numel(xticks)
        if (xticks(k) < 1)
            xTickLabels{k} = sprintf('%2.1f', xticks(k));
        else
            xTickLabels{k} = sprintf('%2.0f', xticks(k));
        end
    end
    
    subplot(1,2,2)
    plot(temporalFreqSupport, ft, 'r-', 'LineWidth', 1.5);
    axis 'square';
    set(gca, 'XScale', 'log', 'XLim', [0.1 100]);
    set(gca, 'YTickLabel', {}, 'XTick', xticks, 'XTickLabel', xTickLabels,'FontSize', 14);
    xlabel('\it temporal frequency (Hz)');
    ylabel('');
    grid on; box on
    pause
    
end

function [coneExcitations, photoCurrents] = retrieveData(stimDescriptor, theContrastLevel, analyzedNoiseInstance, nTrials, eyePosition, resourcesDir)
    fName = coneMosaicResponsesDataFileName(stimDescriptor, theContrastLevel, analyzedNoiseInstance, nTrials, eyePosition, resourcesDir);
    load(fName, 'coneExcitations', 'photoCurrents'); 
end


function energyResponse = computeEnergyResponse(response, nullStimulusResponse, poolingKernelLinear, poolingKernelQuadrature, temporalKernels)
    % Compute modulated mosaic responses by subtracting the MEAN response to the null stimulus
    response = bsxfun(@minus, response, nullStimulusResponse);

    nTrials = size(response,1);
    nTimePoints = size(response,3);
    
    % Convolve with temporal kernel
    plot(squeeze(response(1,1,:)), 'k-'); hold on
    parfor iTrial = 1:nTrials
        tmp = conv2(squeeze(temporalKernels(1,:)), squeeze(response(iTrial,:,:)));
        response(iTrial,:,:) = tmp(:, 1:nTimePoints);
    end
    
    % Compute dot product along all cones for the quadrature pair of pooling kernels
    subunit1Response = squeeze(sum(bsxfun(@times, poolingKernelLinear, response),2));
    subunit2Response = squeeze(sum(bsxfun(@times, poolingKernelQuadrature, response),2));

    % Sum squared outputs
    energyResponse = subunit1Response.^2 + subunit2Response.^2;
    
    
end