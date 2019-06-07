function visualizeEnergyResponses(stimDescriptor, contrastLevels, analyzedNoiseInstance, nTrials, resourcesDir, figNo)

    fname = fullfile(resourcesDir, sprintf('energyResponse_%s_instance_%1.0f_nTrials_%d.mat', stimDescriptor, analyzedNoiseInstance, nTrials));
    load(fname, 'energyConeExcitationResponse', 'energyPhotoCurrentResponse', 'emPathsDegs', 'timeAxis');
   
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 2500 700]);
    
    plotEnergyResponses(stimDescriptor, 'cone excitations', contrastLevels, timeAxis, energyConeExcitationResponse, 0);
    plotEnergyResponses(stimDescriptor, 'photocurrents', contrastLevels, timeAxis, energyPhotoCurrentResponse, 1);
    
end

function plotEnergyResponses(stimDescriptor, signalName, contrastLevels, timeAxis, energyResponse, row)
    nContrasts = size(energyResponse.output,1);
    nTrials = size(energyResponse.output,2);
    for theContrastLevel = 1:nContrasts
        subplot(2, nContrasts, theContrastLevel + row*nContrasts);
        energyResponse.outputMean = squeeze(mean(energyResponse.output,2));
        energyResponse.outputStd = squeeze(std(energyResponse.output,0,2));
        energyResponse.orthoOutputMean = squeeze(mean(energyResponse.orthoOutput,2));
        energyResponse.orthoOutputStd = squeeze(std(energyResponse.orthoOutput,0,2));
        hold on;
        for trialIndex = 1:nTrials
            plot(timeAxis*1000, energyResponse.outputMean(theContrastLevel,:), 'r-', 'LineWidth', 1.5);
            plot(timeAxis*1000, energyResponse.orthoOutputMean(theContrastLevel,:), 'b-', 'LineWidth', 1.5);
             
            plot(timeAxis*1000, energyResponse.outputMean(theContrastLevel,:) + energyResponse.outputStd(theContrastLevel,:), 'r--');
            plot(timeAxis*1000, energyResponse.outputMean(theContrastLevel,:) - energyResponse.outputStd(theContrastLevel,:), 'r--');
            plot(timeAxis*1000, energyResponse.orthoOutputMean(theContrastLevel,:) + energyResponse.orthoOutputStd(theContrastLevel,:), 'b--');
            plot(timeAxis*1000, energyResponse.orthoOutputMean(theContrastLevel,:) - energyResponse.orthoOutputStd(theContrastLevel,:), 'b--');
            
        end
        xlabel('time (msec)');
        ylabel(sprintf('energy response\n(%s)', signalName));
        legend({'standard ori', 'orthogonal ori'});
        title(sprintf('c = %2.3f (%s)', contrastLevels(theContrastLevel),stimDescriptor));
    end
end
    