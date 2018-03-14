function generateFigureForPaper(theFigData, variedParamLegends, variedParamName, fixedParamName, varargin)

    p = inputParser;
    p.addParameter('figureType', 'CSF', @ischar);
    p.addParameter('showBanksPaperIOAcurves', false, @islogical);
    p.addParameter('plotFirstConditionInGray', true, @islogical);
    p.addParameter('inGraphText', '', @ischar);
    p.addParameter('inGraphTextPos', [], @isnumeric);
    
    p.parse(varargin{:});
    
    figureType = p.Results.figureType;
    showBanksPaperIOAcurves = p.Results.showBanksPaperIOAcurves;
    plotFirstConditionInGray = p.Results.plotFirstConditionInGray;
    inGraphText = p.Results.inGraphText;
    inGraphTextPos = p.Results.inGraphTextPos;
    
    for condIndex = 1:numel(theFigData)
        figData = theFigData{condIndex};
        lumIndex = 1;
        referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
        cpd{condIndex} = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
        thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
        contrastSensitivity{condIndex} = 1./(thresholdContrasts*referenceContrast);
    end % condIndex
    
    % Initialize figure
    hFig = figure(1); clf;
    formatFigureForPaper(hFig, 'figureType', figureType);
    
    
    % Displayed colors
    if (plotFirstConditionInGray)
        colors(1,:) = [0.5 0.5 0.5];
        if (numel(theFigData)>1)
            colors(2:numel(theFigData),:) = brewermap(numel(theFigData)-1, 'Set1');
        end
    else
        colors = brewermap(numel(theFigData), 'Set1');
    end
    
    faceColor = [0.3 0.3 0.3];
    
    hold on;

    if (showBanksPaperIOAcurves)
        banksFactor = 1;
        plot(figData.D(:,1),figData.D(:,2)*banksFactor,'-','Color', squeeze(colors(1,:)), 'LineWidth',2);
        plot(figData.C(:,1),figData.C(:,2)*banksFactor,'-','Color', squeeze(colors(2,:)), 'LineWidth',2);
        plot(figData.E(:,1),figData.E(:,2)*banksFactor,'-','Color', squeeze(colors(3,:)), 'LineWidth',2);
        for condIndex = 1:numel(theFigData)
            plot(cpd{condIndex}, contrastSensitivity{condIndex}, 'o', 'Color', squeeze(colors(condIndex,:)), ...
            'MarkerEdgeColor', squeeze(colors(condIndex,:)), ...
            'MarkerFaceColor', min(1, squeeze(colors(condIndex,:)) + faceColor), ...
            'MarkerSize',12,'LineWidth',2);
        end
    else
        for condIndex = 1:numel(theFigData)
            plot(cpd{condIndex}, contrastSensitivity{condIndex}, 'o-', 'Color', squeeze(colors(condIndex,:)), ...
            'MarkerEdgeColor', squeeze(colors(condIndex,:)), ...
            'MarkerFaceColor', min(1, squeeze(colors(condIndex,:)) + faceColor), ...
            'MarkerSize',12,'LineWidth',2);
        end
    end
    
    % Add legend
    hL = legend(variedParamLegends);
    
    % Add text
    if (~isempty(inGraphText))
        if (isempty(inGraphTextPos))
            inGraphTextPos(1) = 1.1;
            inGraphTextPos(2) = 9000;
        end
        t = text(inGraphTextPos(1), inGraphTextPos(2), inGraphText);
    end

    
    % Format figure
    formatFigureForPaper(hFig, 'theAxes', gca, 'theLegend', hL, 'theText', t, 'figureType', 'CSF');
    
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    figureName = fullfile(exportsDir, sprintf('%sVary%s.pdf', variedParamName, fixedParamName));
    NicePlot.exportFigToPDF(figureName, hFig, 300);
end

