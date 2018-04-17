function generateFigureForPaper(theFigData, variedParamLegends, variedParamName, fixedParamName, varargin)

    p = inputParser;
    p.addParameter('figureType', 'CSF', @ischar);
    p.addParameter('showBanksPaperIOAcurves', false, @islogical);
    p.addParameter('plotFirstConditionInGray', true, @islogical);
    p.addParameter('plotRatiosOfOtherConditionsToFirst', false, @islogical);
    p.addParameter('theRatioLims', []);
    p.addParameter('theRatioTicks', []);
    p.addParameter('inGraphText', '', @ischar);
    p.addParameter('inGraphTextPos', [], @isnumeric);
    p.addParameter('inGraphTextFontSize', [], @isnumeric);
    p.parse(varargin{:});
    
    figureType = p.Results.figureType;
    showBanksPaperIOAcurves = p.Results.showBanksPaperIOAcurves;
    plotFirstConditionInGray = p.Results.plotFirstConditionInGray;
    inGraphText = p.Results.inGraphText;
    inGraphTextPos = p.Results.inGraphTextPos;
    inGraphTextFontSize = p.Results.inGraphTextFontSize;
    plotRatiosOfOtherConditionsToFirst = p.Results.plotRatiosOfOtherConditionsToFirst;
    theRatioLims = p.Results.theRatioLims;
    theRatioTicks = p.Results.theRatioTicks;
    
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
    [theAxes, theRatioAxes] = formatFigureForPaper(hFig, ...
        'figureType', figureType, ...
        'plotRatiosOfOtherConditionsToFirst', plotRatiosOfOtherConditionsToFirst);
    
    
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
    
    hold(theAxes(1), 'on');

    if (showBanksPaperIOAcurves)
        %banksFactor = 1;
        %plot(theAxes,figData.D(:,1),figData.D(:,2)*banksFactor,'-','Color', squeeze(colors(1,:)), 'LineWidth',2);
        %plot(theAxes,figData.C(:,1),figData.C(:,2)*banksFactor,'-','Color', squeeze(colors(2,:)), 'LineWidth',2);
        %plot(theAxes,figData.E(:,1),figData.E(:,2)*banksFactor,'-','Color', squeeze(colors(3,:)), 'LineWidth',2);
        
        plot(theAxes,figData.BanksCSF('34 cd/m2').x,figData.BanksCSF('34 cd/m2').y,'-','Color', squeeze(colors(1,:)), 'LineWidth',2);
        plot(theAxes,figData.BanksCSF('340 cd/m2').x,figData.BanksCSF('340 cd/m2').y,'-','Color', squeeze(colors(2,:)), 'LineWidth',2);
        plot(theAxes,figData.BanksCSF('3.4 cd/m2').x,figData.BanksCSF('3.4 cd/m2').y,'-','Color', squeeze(colors(3,:)), 'LineWidth',2);
        
        for condIndex = 1:numel(theFigData)
            plot(theAxes, cpd{condIndex}, contrastSensitivity{condIndex}, 'o', 'Color', squeeze(colors(condIndex,:)), ...
            'MarkerEdgeColor', squeeze(colors(condIndex,:)), ...
            'MarkerFaceColor', min(1, squeeze(colors(condIndex,:)) + faceColor), ...
            'MarkerSize',12,'LineWidth',2);
        end
    else
        for condIndex = 1:numel(theFigData)
            plot(theAxes, cpd{condIndex}, contrastSensitivity{condIndex}, 'o-', 'Color', squeeze(colors(condIndex,:)), ...
            'MarkerEdgeColor', max(0, squeeze(colors(condIndex,:))-faceColor), ...
            'MarkerFaceColor', min(1, squeeze(colors(condIndex,:)) + faceColor), ...
            'MarkerSize',12,'LineWidth',2);
        end
    end
    
    % Add legend
    hL = legend(theAxes, variedParamLegends);
    
    % Add text
    if (~isempty(inGraphText))
        if (isempty(inGraphTextPos))
            if (plotRatiosOfOtherConditionsToFirst)
                inGraphTextPos(1) = 1.1;
                inGraphTextPos(2) = 8500;
            else
                inGraphTextPos(1) = 1.1;
                inGraphTextPos(2) = 9000;
            end
        end
        t = text(theAxes, inGraphTextPos(1), inGraphTextPos(2), inGraphText);
    else
        t = '';
    end
    hold(theAxes, 'off');
    
    if (plotRatiosOfOtherConditionsToFirst)
        
        hold(theRatioAxes, 'on');
        refCPD = cpd{1};
        refSensitivity = contrastSensitivity{1};
        for condIndex = 2:numel(theFigData)
            condCPD = cpd{condIndex};
            condSensitivity = contrastSensitivity{condIndex};
            ratios = zeros(1,numel(refCPD));
            for k = 1:numel(condCPD)
                idx = find(condCPD(k) == refCPD);
                if (~isempty(idx))
                    ratios(k) = condSensitivity(idx)./refSensitivity(idx);
                end
            end
            plot(theRatioAxes, refCPD, ratios, ...
                'ko-', 'Color', squeeze(colors(condIndex,:)), ...
                'MarkerEdgeColor', squeeze(colors(condIndex,:)), ...
                'MarkerFaceColor', min(1, squeeze(colors(condIndex,:)) + faceColor), ...
                'MarkerSize',12,'LineWidth',2);
        end
        hold(theRatioAxes, 'off');
    end
    
    
    % Format figure
    formatFigureForPaper(hFig, ...
        'figureType', 'CSF', ...
        'plotRatiosOfOtherConditionsToFirst', plotRatiosOfOtherConditionsToFirst, ...
        'theAxes', theAxes, ...
        'theRatioAxes', theRatioAxes, ...
        'theRatioLims', theRatioLims, ...
        'theRatioTicks', theRatioTicks, ...
        'theLegend', hL, ...
        'theText', t, ...
        'theTextFontSize', inGraphTextFontSize);
    
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    figureName = fullfile(exportsDir, sprintf('%sVary%s.pdf', variedParamName, fixedParamName));
    NicePlot.exportFigToPDF(figureName, hFig, 300);
end

