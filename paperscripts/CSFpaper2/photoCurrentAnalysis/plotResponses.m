function plotResponses(modelResponse, modelResponseOpositePolarity, adaptationLevel, contrastLevel, ...
    pulseDurationSeconds, photoCurrentRange, showSNR, showSNRcomponents, coneExcitationRange, figNo)
 
    if (showSNR)
        [coneExcitation.SNR, ...
         coneExcitation.noiseEstimationLatency, ...
         coneExcitation.peakEstimationLatency, ...
         coneExcitation.modulationPeak, ...
         coneExcitation.noiseSigma, ...
         coneExcitation.responseAtPeakSigma, ...
         coneExcitation.responseAtPeak] = ...
            computeSNR(modelResponse.noisyConeExcitationTimeAxis, modelResponse.noisyConeExcitations, modelResponse.meanConeExcitationCountSignal);
        
        
        [photoCurrent.SNR, ...
         photoCurrent.noiseEstimationLatency, ...
         photoCurrent.peakEstimationLatency, ...
         photoCurrent.modulationPeak, ...
         photoCurrent.noiseSigma, ...
         photoCurrent.responseAtPeakSigma, ...
         photoCurrent.responseAtPeak] = ...
            computeSNR(modelResponse.timeAxis, modelResponse.noisyMembraneCurrents, modelResponse.membraneCurrent);
        
        if (~isempty(modelResponseOpositePolarity))
            [coneExcitationOpositePolarity.SNR, ...
             coneExcitationOpositePolarity.noiseEstimationLatency, ....
             coneExcitationOpositePolarity.peakEstimationLatency, ...
             coneExcitationOpositePolarity.modulationPeak, ...
             coneExcitationOpositePolarity.noiseSigma, ...
             coneExcitationOpositePolarity.responseAtPeakSigma, ...
             coneExcitationOpositePolarity.responseAtPeak] = ...
                computeSNR(modelResponseOpositePolarity.noisyConeExcitationTimeAxis, modelResponseOpositePolarity.noisyConeExcitations, modelResponseOpositePolarity.meanConeExcitationCountSignal);
            
            [photoCurrentOpositePolarity.SNR, ...
             photoCurrentOpositePolarity.noiseEstimationLatency, ...
             photoCurrentOpositePolarity.peakEstimationLatency, ...
             photoCurrentOpositePolarity.modulationPeak, ...
             photoCurrentOpositePolarity.noiseSigma, ...
             photoCurrentOpositePolarity.responseAtPeakSigma, ...
             photoCurrentOpositePolarity.responseAtPeak] = ...
                computeSNR(modelResponseOpositePolarity.timeAxis, modelResponseOpositePolarity.noisyMembraneCurrents, modelResponseOpositePolarity.membraneCurrent);
        
        end
        
    end
    
    if (isempty(coneExcitationRange))
        adaptingExcitationCount = pulseDurationSeconds * adaptationLevel;
        coneExcitationRange = [0 2*adaptingExcitationCount];
        [coneExcitationTicks, coneExcitationTickLabels] = coneExcTicksAndLabels(2*adaptingExcitationCount);
        
    else
        [coneExcitationTicks, coneExcitationTickLabels] = coneExcTicksAndLabels(coneExcitationRange(2));
    end
    
    
    if (photoCurrentRange(2)-photoCurrentRange(1)>=50)
        photoCurrentTicks = [-90:10:0];
    elseif (photoCurrentRange(2)-photoCurrentRange(1)>=20)
        photoCurrentTicks = [-90:5:0];
    elseif (photoCurrentRange(2)-photoCurrentRange(1)>=10)
        photoCurrentTicks = [-90:2:0];
    elseif (photoCurrentRange(2)-photoCurrentRange(1)>=5)
        photoCurrentTicks = [-90:1:0];
    elseif (photoCurrentRange(2)-photoCurrentRange(1)>=1)
        photoCurrentTicks = [-90:0.5:0];
    elseif (photoCurrentRange(2)-photoCurrentRange(1)>=0.5)
        photoCurrentTicks = [-90:0.1:0];
    else
        photoCurrentTicks = [-90:0.05:0];
    end
    
    timeLims = [-0.1 0.8];
    timeTicks = [0:0.2:1.0];
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 350 460], 'Color', [1 1 1]);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', 1, ...
       'rowsNum', 2, ...
       'heightMargin',   0.06, ...
       'widthMargin',    0.03, ...
       'leftMargin',     0.17, ...
       'rightMargin',    0.02, ...
       'bottomMargin',   0.14, ...
       'topMargin',      0.04);
    
    extraTime = 0.5;
    dt = modelResponse.noisyConeExcitationTimeAxis(2)-modelResponse.noisyConeExcitationTimeAxis(1);
    idx = find(modelResponse.noisyConeExcitationTimeAxis>=modelResponse.noisyConeExcitationTimeAxis(end)-extraTime);
    negTime = modelResponse.noisyConeExcitationTimeAxis(idx)-modelResponse.noisyConeExcitationTimeAxis(end)-dt;
    
    % Add response component before time 0
    modelResponse.noisyConeExcitationTimeAxis = cat(2,negTime, modelResponse.noisyConeExcitationTimeAxis);
    modelResponse.noisyConeExcitations = cat(2,modelResponse.noisyConeExcitations(:,idx), modelResponse.noisyConeExcitations);
    modelResponse.meanConeExcitationCountSignal = cat(2, modelResponse.meanConeExcitationCountSignal(:,idx), modelResponse.meanConeExcitationCountSignal);
    modelResponse.noisyConeExcitationTimeAxis = modelResponse.noisyConeExcitationTimeAxis + dt;
    
    if (~isempty(modelResponseOpositePolarity))
        modelResponseOpositePolarity.noisyConeExcitationTimeAxis = cat(2,negTime, modelResponseOpositePolarity.noisyConeExcitationTimeAxis);
        modelResponseOpositePolarity.noisyConeExcitations = cat(2,modelResponseOpositePolarity.noisyConeExcitations(:,idx), modelResponseOpositePolarity.noisyConeExcitations);
        modelResponseOpositePolarity.meanConeExcitationCountSignal = cat(2, modelResponseOpositePolarity.meanConeExcitationCountSignal(:,idx), modelResponseOpositePolarity.meanConeExcitationCountSignal);
        modelResponseOpositePolarity.noisyConeExcitationTimeAxis = modelResponseOpositePolarity.noisyConeExcitationTimeAxis + dt;
    end
    
    idx = find(modelResponse.timeAxis>modelResponse.timeAxis(end)-extraTime);
    dt = modelResponse.timeAxis(2)-modelResponse.timeAxis(1);
    negTime = modelResponse.timeAxis(idx) - modelResponse.timeAxis(end)-dt;
    
    % Add response component before time0
    modelResponse.timeAxis(modelResponse.timeAxis<0) = nan;
    modelResponse.timeAxis = cat(2, negTime, modelResponse.timeAxis);
    modelResponse.noisyMembraneCurrents = cat(2, modelResponse.noisyMembraneCurrents(:,idx), modelResponse.noisyMembraneCurrents);
    modelResponse.membraneCurrent = cat(2, modelResponse.membraneCurrent(idx), modelResponse.membraneCurrent);
    if (~isempty(modelResponseOpositePolarity))
        modelResponseOpositePolarity.timeAxis(modelResponseOpositePolarity.timeAxis<0) = nan;
        modelResponseOpositePolarity.timeAxis = cat(2, negTime, modelResponseOpositePolarity.timeAxis);
        modelResponseOpositePolarity.noisyMembraneCurrents = cat(2, modelResponseOpositePolarity.noisyMembraneCurrents(:,idx), modelResponseOpositePolarity.noisyMembraneCurrents);
        modelResponseOpositePolarity.membraneCurrent = cat(2, modelResponseOpositePolarity.membraneCurrent(idx), modelResponseOpositePolarity.membraneCurrent);
    end
    
    % Display up to 500 instances
    if (isempty(modelResponseOpositePolarity))
        displayedInstancesNum = 500;
    else
        % 250 instances from one polarity, 250 from the other
        displayedInstancesNum = 250;
    end
    displayedInstancesNum = min([displayedInstancesNum  size(modelResponse.noisyConeExcitations,1)]);
    
    modelResponse.noisyConeExcitations = modelResponse.noisyConeExcitations(1:displayedInstancesNum,:);
    modelResponse.noisyMembraneCurrents = modelResponse.noisyMembraneCurrents(1:displayedInstancesNum,:);
    if (~isempty(modelResponseOpositePolarity))
        modelResponseOpositePolarity.noisyConeExcitations = modelResponseOpositePolarity.noisyConeExcitations(1:displayedInstancesNum,:);
        modelResponseOpositePolarity.noisyMembraneCurrents = modelResponseOpositePolarity.noisyMembraneCurrents(1:displayedInstancesNum,:);
    end

    
    subplot('Position', subplotPosVectors(1,1).v);
    plotResponseInstances = true;
    theYLabel = sprintf('\\it cone excitations (R*/c/%0.2fs)', pulseDurationSeconds);
    stairsComboPlot('g', modelResponse.noisyConeExcitationTimeAxis, modelResponse.noisyConeExcitations, modelResponse.meanConeExcitationCountSignal, ...
        coneExcitationRange, coneExcitationTicks, coneExcitationTickLabels, timeLims, timeTicks, false, true, false, theYLabel, plotResponseInstances);
    
    if (~isempty(modelResponseOpositePolarity))
        stairsComboPlot('r', modelResponseOpositePolarity.noisyConeExcitationTimeAxis, modelResponseOpositePolarity.noisyConeExcitations, ...
            modelResponseOpositePolarity.meanConeExcitationCountSignal, coneExcitationRange, coneExcitationTicks, coneExcitationTickLabels, timeLims, timeTicks, false, false, false, theYLabel, plotResponseInstances);
    end
    
    plotResponseInstances = false;
    stairsComboPlot('g', modelResponse.noisyConeExcitationTimeAxis, modelResponse.noisyConeExcitations, modelResponse.meanConeExcitationCountSignal, ...
        coneExcitationRange, coneExcitationTicks, coneExcitationTickLabels, timeLims, timeTicks, false, true, false, theYLabel, plotResponseInstances);
    if (~isempty(modelResponseOpositePolarity))
        stairsComboPlot('r', modelResponseOpositePolarity.noisyConeExcitationTimeAxis, modelResponseOpositePolarity.noisyConeExcitations, ...
            modelResponseOpositePolarity.meanConeExcitationCountSignal, coneExcitationRange, coneExcitationTicks, coneExcitationTickLabels, timeLims, timeTicks, false, false, false, theYLabel, plotResponseInstances);
    end
    
    % Show SNR components
    if (showSNR)
        if (showSNRcomponents)
            plot([0 0], modelResponse.meanConeExcitationCountSignal(end) + [0 -coneExcitation.modulationPeak], 'g-', 'LineWidth', 1.5);
            plot(coneExcitation.peakEstimationLatency+pulseDurationSeconds/2*[1 1], coneExcitation.responseAtPeak + coneExcitation.responseAtPeakSigma*[-1 1], 'gs-', 'LineWidth', 1.5); 
            plot(-0.05*[1 1], modelResponse.meanConeExcitationCountSignal(end) + coneExcitation.noiseSigma*[-1 1], 'y-', 'LineWidth', 1.5); 

            if (~isempty(modelResponseOpositePolarity))
                plot([0 0], modelResponseOpositePolarity.meanConeExcitationCountSignal(end) + [0 coneExcitationOpositePolarity.modulationPeak], 'r-', 'LineWidth', 1.5);
                plot(coneExcitationOpositePolarity.peakEstimationLatency+pulseDurationSeconds/2*[1 1], coneExcitationOpositePolarity.responseAtPeak + coneExcitationOpositePolarity.responseAtPeakSigma*[-1 1], 'rs-', 'LineWidth', 1.5); 
            end
        end
        
        if (~isempty(modelResponseOpositePolarity))
            title(sprintf('SNR(dec/inc): %2.2f/%2.2f', coneExcitation.SNR, coneExcitationOpositePolarity.SNR), 'FontWeight', 'normal');
        else
            title(sprintf('SNR: %2.2f/%2.2f', coneExcitation.SNR), 'FontWeight', 'normal');
        end
    end
    
    subplot('Position', subplotPosVectors(2,1).v);
    plotResponseInstances = true;
    linesComboPlot('g', modelResponse.timeAxis, modelResponse.noisyMembraneCurrents', modelResponse.membraneCurrent, ...
        photoCurrentRange, photoCurrentTicks, timeLims, timeTicks, false, true, true, plotResponseInstances);
    if (~isempty(modelResponseOpositePolarity))
        linesComboPlot('r', modelResponseOpositePolarity.timeAxis, modelResponseOpositePolarity.noisyMembraneCurrents', ...
            modelResponseOpositePolarity.membraneCurrent, photoCurrentRange, photoCurrentTicks, timeLims, timeTicks, false, true, true, plotResponseInstances);
    end
    
    plotResponseInstances = false;
    linesComboPlot('g', modelResponse.timeAxis, modelResponse.noisyMembraneCurrents', modelResponse.membraneCurrent, ...
        photoCurrentRange, photoCurrentTicks, timeLims, timeTicks, true, true, true, plotResponseInstances);
    if (~isempty(modelResponseOpositePolarity))
        linesComboPlot('r', modelResponseOpositePolarity.timeAxis, modelResponseOpositePolarity.noisyMembraneCurrents', ...
            modelResponseOpositePolarity.membraneCurrent, photoCurrentRange, photoCurrentTicks, timeLims, timeTicks, true, true, true, plotResponseInstances);
    end
    
    % Show SNR components
    if (showSNR)
        if (showSNRcomponents)
            plot([0 0], modelResponse.membraneCurrent(end) + [0 -photoCurrent.modulationPeak], 'g-', 'LineWidth', 1.5);
            plot(photoCurrent.peakEstimationLatency+[0 0], photoCurrent.responseAtPeak + photoCurrent.responseAtPeakSigma*[-1 1], 'gs-', 'LineWidth', 1.5); 
            plot(-0.05*[1 1], modelResponse.membraneCurrent(end) + photoCurrent.noiseSigma*[-1 1], 'y-', 'LineWidth', 1.5); 

            if (~isempty(modelResponseOpositePolarity))
                plot([0 0], modelResponseOpositePolarity.membraneCurrent(end) + [0 photoCurrentOpositePolarity.modulationPeak], 'r-', 'LineWidth', 1.5);
                plot(photoCurrentOpositePolarity.peakEstimationLatency+[0 0], photoCurrentOpositePolarity.responseAtPeak + photoCurrentOpositePolarity.responseAtPeakSigma*[-1 1], 'rs-', 'LineWidth', 1.5); 
            end
        end
    
        if (~isempty(modelResponseOpositePolarity))
            title(sprintf('SNR(dec/inc): %2.2f/%2.2f', photoCurrent.SNR, photoCurrentOpositePolarity.SNR), 'FontWeight', 'normal');
        else
            title(sprintf('SNR: %2.2f/%2.2f', photoCurrent.SNR), 'FontWeight', 'normal');
        end
    end
    
    drawnow;
    if (isempty(modelResponseOpositePolarity))
        pdfFileName = sprintf('Responses_adapt_%2.0f_contrast_%2.3f_%2.0fmsec.pdf', adaptationLevel, contrastLevel, pulseDurationSeconds*1000);
    else
        pdfFileName = sprintf('Responses_adapt_%2.0f_contrast_pm%2.3f_%2.0fmsec.pdf', adaptationLevel, abs(contrastLevel), pulseDurationSeconds*1000);
    end
    
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300)
    
end

function stairsComboPlot(color, timeAxis, responseInstances, meanResponse, ...
    responseLims, responseTicks, responseTickLabels, timeLims, timeTicks, ...
    showXLabel, showYLabel, showXTickLabel, theYLabel, plotResponseInstances)
%
    if (plotResponseInstances)
        [stairsTime, stairsResponse] = stairsCoords(timeAxis, responseInstances);
        hPlot = plot(stairsTime, stairsResponse,  '-', 'Color', 'k'); hold on;
        for k = 1:numel(hPlot)
            hPlot(k).Color(4) = 0.3;  % 5% transparent
        end
        hold on;
    end
    
    [stairsTime, stairsResponse] = stairsCoords(timeAxis, meanResponse);
    plot(stairsTime, stairsResponse, '-', 'Color', color,  'LineWidth', 2.0);
    
    if (showYLabel)
        ylabel(theYLabel);
    end
    if (showXLabel)
        xlabel('\it time (sec)');
    end
    
    grid on
    set(gca, 'FontSize', 14);
    set(gca, 'XLim', timeLims, 'XTick', timeTicks, 'YLim', responseLims, ...
        'YTick', responseTicks, 'YTickLabel', responseTickLabels);
    
    if (~showXTickLabel)
        set(gca, 'XTickLabel', {});
    end
end

function linesComboPlot(color, timeAxis, responseInstances, meanResponse, responseRange, responseTicks, timeLims, timeTicks, showXLabel, showYLabel, showXTickLabel, plotResponseInstances)
    if (plotResponseInstances)
        hPlot = plot(timeAxis, responseInstances, '-', 'Color', 'k');
        for k = 1:numel(hPlot)
            hPlot(k).Color(4) = 0.05;  % 5% transparent
        end
        hold on;
    end
    
    plot(timeAxis, meanResponse, '-', 'Color', color,  'LineWidth', 2.0);
    grid on;
    set(gca,  'XLim', timeLims, 'XTick', timeTicks, 'YLim', responseRange, 'YTick', responseTicks);
    if (showYLabel)
        ylabel('\it photocurrent (pA)');  
    end
    if (showXLabel)
        xlabel('\it time (sec)');
    end
    if (~showXTickLabel)
        set(gca, 'XTickLabel', {});
    end
    
    set(gca, 'FontSize', 14);
end


function [xStairs, yStairs] = stairsCoords(x,y)
    xStairs(1,1) = x(1);
    yStairs(:,1) = y(:,1);
    for k = 1:numel(x)-1
        xStairs = cat(2,xStairs, x(1,k));
        yStairs = cat(2,yStairs, y(:,k));
        
        xStairs = cat(2,xStairs, x(1,k));
        yStairs = cat(2,yStairs, y(:,k+1));
    end
    xStairs = cat(2, xStairs, x(1,k));
    yStairs = cat(2, yStairs, y(:,k));
end

function [coneExcitationTicks, coneExcitationTickLabels] = coneExcTicksAndLabels(adaptingExcitationCount)
    if (adaptingExcitationCount <=50)
        coneExcitationTicks = 0:10:100;
        coneExcitationTickLabels = {'0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'};
    elseif (adaptingExcitationCount <=100) 
        coneExcitationTicks = 0:25:200;
        coneExcitationTickLabels = {'0', '25', '50', '70', '100', '125', '150', '170', '200'};
    elseif (adaptingExcitationCount <=200) 
        coneExcitationTicks = 0:50:400;
        coneExcitationTickLabels = {'0', '50', '100', '150', '200', '250', '300', '350', '400'};
    elseif (adaptingExcitationCount <=400)
        coneExcitationTicks = 0:100:700;
        coneExcitationTickLabels = {'0', '0.1k', '0.2k', '0.3k', '0.4k', '0.5k', '0.6k', '0.7k'};
    elseif (adaptingExcitationCount <=800)
        coneExcitationTicks = 0:200:1500;
        coneExcitationTickLabels = {'0', '0.2k', '0.4k', '0.6k', '0.8k', '1.0k', '1.2k', '1.4k'};
    elseif (adaptingExcitationCount <=1600)
        coneExcitationTicks = 0:500:3000;
        coneExcitationTickLabels = {'0', '0.5k', '1.0k', '1.5k', '2.0k', '2.5k', '3.0k'};
    elseif (adaptingExcitationCount <=3200)
        coneExcitationTicks = 0:500:6000;
        coneExcitationTickLabels = {'0', '0.5k', '1.0k', '1.5k', '2.0k', '2.5k', '3.0k', '3.5k', '4.0k', '4.5k', '5.0k', '5.5k', '6.0k'};
    elseif (adaptingExcitationCount <=6400)
        coneExcitationTicks = 0:1000:15000;
        coneExcitationTickLabels = {'0', '1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '10k', '10k', '11k'};
    else
        coneExcitationTicks = 0:5000:60000;
        coneExcitationTickLabels = {'0', '5k', '10k', '15k', '20k', '25k', '30k', '35k', '40k', '45k', '50k', '55k', '60k'};            
    end
end

