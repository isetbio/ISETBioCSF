function visualizeEnergyResponses(stimDescriptor, signalName, ...
    standardOriStimulusResponse, standardOriStimulusOrthoResponse, ...
    orthogonalOriStimulusResponse, orthogonalOriStimulusOrthoResponse, ...
    contrastLevels, timeAxis, figNo)

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [61 22 2500 700]);
    signalName
    %plotEnergyResponses(stimDescriptor, 'stardard orientation', signalName, contrastLevels, timeAxis, ...
    %    standardOriStimulusResponse, standardOriStimulusOrthoResponse, 0);
    %plotEnergyResponses(stimDescriptor, 'orthogonal orientation', signalName, contrastLevels, timeAxis, ...
    %    orthogonalOriStimulusResponse, orthogonalOriStimulusOrthoResponse, 1);
    plotEnergyResponseCombo(standardOriStimulusResponse, standardOriStimulusOrthoResponse, ...
        orthogonalOriStimulusResponse, orthogonalOriStimulusOrthoResponse, signalName, contrastLevels, stimDescriptor, 2);
    pause
end

function plotEnergyResponseCombo(energyResponseStandardOriStimulusOutput, energyResponseStandardOriStimulusOrthoOutput, ...
    energyResponseOrthogonalOriStimulusOutput, energyResponseOrthogonalOriStimulusOrthoOutput, signalName, contrastLevels, stimDescriptor, row)
    nContrasts = size(energyResponseStandardOriStimulusOutput,1);
    nTrials = size(energyResponseStandardOriStimulusOutput,2);
    nTimeBins = size(energyResponseStandardOriStimulusOutput,3);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', nContrasts, ...
       'heightMargin',  0.01, ...
       'widthMargin',   0.01, ...
       'leftMargin',    0.03, ...
       'rightMargin',   0.01, ...
       'bottomMargin',  0.05, ...
       'topMargin',     0.03);
   
    maxAll = 0.8*max([...
        max(energyResponseStandardOriStimulusOutput(:))
        max(energyResponseStandardOriStimulusOrthoOutput(:))
        max(energyResponseOrthogonalOriStimulusOutput(:))
        max(energyResponseOrthogonalOriStimulusOrthoOutput(:))]);
    
    for theContrastLevel = 1:nContrasts
        subplot('Position', subplotPosVectors(1, theContrastLevel).v);
        
        r1 = energyResponseStandardOriStimulusOutput(theContrastLevel,1:nTrials,1:nTimeBins);
        r2 = energyResponseStandardOriStimulusOrthoOutput(theContrastLevel,1:nTrials,1:nTimeBins);
        r3 = energyResponseOrthogonalOriStimulusOutput(theContrastLevel,1:nTrials,1:nTimeBins);
        r4 = energyResponseOrthogonalOriStimulusOrthoOutput(theContrastLevel,1:nTrials,1:nTimeBins);
        
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
        xlabel('\it 0 deg mechanism energy');
        axis 'square'
        set(gca, 'FontSize', 16, 'XLim', [0 maxAll], 'XTick', 0:0.01:0.1, 'YTick', 0:0.01:0.1, 'YLim', [0 maxAll]);
        if (theContrastLevel>1)
            set(gca, 'YTickLabel', {});
        else
            ylabel('\it 90 deg mechanism energy'); 
        end
        
        legend({'0 deg stimulus', '90 deg stimulus'});
        axis 'square';
        title(sprintf('c = %2.2f%% (%s,%s)', contrastLevels(theContrastLevel)*100, stimDescriptor, signalName));
        drawnow
    end
end

function plotEnergyResponses(stimDescriptor, stimOrientation, signalName, contrastLevels, timeAxis, energyResponseOutput, energyResponseOrthoOutput, row)
    nContrasts = size(energyResponseOutput,1);
    nTrials = size(energyResponseOutput,2);
    for theContrastLevel = 1:nContrasts
        subplot(3, nContrasts, theContrastLevel + row*nContrasts);
        energyResponse.outputMean = squeeze(mean(energyResponseOutput,2));
        energyResponse.outputStd = squeeze(std(energyResponseOutput,0,2));
        energyResponse.orthoOutputMean = squeeze(mean(energyResponseOrthoOutput,2));
        energyResponse.orthoOutputStd = squeeze(std(energyResponseOrthoOutput,0,2));
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
        drawnow;
    end
end
    