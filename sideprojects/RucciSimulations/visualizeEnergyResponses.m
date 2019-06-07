function visualizeEnergyResponses(stimDescriptor, signalName, responseStandardOriStimulus, responseOrthogonalOriStimulus, contrastLevels, timeAxis, figNo)

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 2500 700]);
    
    plotEnergyResponses(stimDescriptor, 'stardard orientation', signalName, contrastLevels, timeAxis, responseStandardOriStimulus, 0);
    plotEnergyResponses(stimDescriptor, 'orthogonal orientation', signalName, contrastLevels, timeAxis, responseOrthogonalOriStimulus, 1);
    plotEnergyResponseCombo(responseStandardOriStimulus, responseOrthogonalOriStimulus, signalName, 2);
end

function plotEnergyResponseCombo(energyResponseStandardOriStimulus, energyResponseOrthogonalOriStimulus, signalName, row)
    nContrasts = size(energyResponseStandardOriStimulus.output,1);
    nTrials = size(energyResponseStandardOriStimulus.output,2);
    nTimeBins = size(energyResponseStandardOriStimulus.output,3);
    
    for theContrastLevel = 1:nContrasts
        subplot(3, nContrasts, theContrastLevel + row*nContrasts);
        
        r1 = energyResponseStandardOriStimulus.output(theContrastLevel,1:nTrials,1:nTimeBins);
        r2 = energyResponseStandardOriStimulus.orthoOutput(theContrastLevel,1:nTrials,1:nTimeBins);
        r3 = energyResponseOrthogonalOriStimulus.output(theContrastLevel,1:nTrials,1:nTimeBins);
        r4 = energyResponseOrthogonalOriStimulus.orthoOutput(theContrastLevel,1:nTrials,1:nTimeBins);
        
        r1 = r1(:);
        r2 = r2(:);
        r3 = r3(:);
        r4 = r4(:);
        
        rr = [r1; r2; r3; r4];
        responseRange = prctile(rr, [1 99]);
        
        plot(r1,r2,'r.'); hold on;    
        plot(r3, r4, 'b.');
        plot(responseRange(1)*[1 1], responseRange(2)*[1 1], 'k-');
        
        set(gca, 'XLim', responseRange, 'YLim', responseRange);
        xlabel('standard ori mechanism');
        ylabel('orthogonal ori mechanism'); 
        legend({'standard ori stimulus', 'orthogonal ori stimulus'});
        axis 'square';
        title(signalName);
    end
end

function plotEnergyResponses(stimDescriptor, stimOrientation, signalName, contrastLevels, timeAxis, energyResponse, row)
    nContrasts = size(energyResponse.output,1);
    nTrials = size(energyResponse.output,2);
    for theContrastLevel = 1:nContrasts
        subplot(3, nContrasts, theContrastLevel + row*nContrasts);
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
        legend({'standard ori mechanism', 'orthogonal ori mechanism'});
        title(sprintf('c = %2.3f (%s, %s)', contrastLevels(theContrastLevel),stimDescriptor, stimOrientation));
    end
end
    