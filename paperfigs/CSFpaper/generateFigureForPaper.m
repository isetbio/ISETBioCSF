function generateFigureForPaper(theFigData, variedParamLegends, variedParamName, fixedParamName, varargin)

    p = inputParser;
    p.addParameter('figureType', 'CSF', @ischar);
    p.parse(varargin{:});
    
    figureType = p.Results.figureType;
    
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
    colors(1,:) = [0.5 0.5 0.5];
    if (numel(theFigData)>1)
        colors(2:numel(theFigData),:) = brewermap(numel(theFigData)-1, 'Set1');
    end
    faceColor = [0.3 0.3 0.3];
    
    hold on;
    for condIndex = 1:numel(theFigData)
        plot(cpd{condIndex}, contrastSensitivity{condIndex}, 'o-', 'Color', squeeze(colors(condIndex,:)), ...
        'MarkerEdgeColor', squeeze(colors(condIndex,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(condIndex,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    end
    
    % Add legend
    hL = legend(variedParamLegends);
    
    % Format figure
    formatFigureForPaper(hFig, 'theAxes', gca, 'theLegend', hL, 'figureType', 'CSF');
    
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    figureName = fullfile(exportsDir, sprintf('%sVary%s.pdf', variedParamName, fixedParamName));
    NicePlot.exportFigToPDF(figureName, hFig, 300);
end

