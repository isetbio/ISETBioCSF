function generateCrowellBanksFigure(data)
    
    performanceClassifier = 'svmV1FilterEnsemble';
    performanceClassifier = 'svmV1FilterBank';
    
    if notDefined('data')
        dataFilename = sprintf('data_mlpt.mat', strrep(performanceClassifier, ' ', '_'));
        load(fullfile('/Users/nicolas/Desktop',dataFilename), 'data');
        referenceData = data;
        
        dataFilename = sprintf('data_%s.mat', strrep(performanceClassifier, ' ', '_'));
        load(fullfile('/Users/nicolas/Desktop',dataFilename), 'data');
    end
    
        
    
    csTicks = [2 5 10 20 50 100 200 500 1000 2000 5000 10000];
    csLims = [1.5 6000];
    sfTicks = [1 2 5 10 20 50 100];
    sfLims  = [1.5 80];
    dx1 = 0.2; dx2 = 10; dy1 = 0.5; dy2 = 1000;
    theRatioTicks = [0.01 0.03 0.1 0.3 1];
    theRatioLims = [0.01 1.0];
    
    hFig = figure(1234); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 800 340]);
    
    
    
    markerSize = 6;
    lineWidth = 1.15;
    fontSize = 9;
    legendFontSize = 7;
    textFontSize = 12;
    labelFontSize = 10;
    axesLineWidth = 0.25;
           
    load('CrowellBanksSubjects.mat');
    
    % Plot the computed data

    for patchSizeIndex  = 1:numel(data.patchSize2SigmaCycles)
        
        ratioAxes = subplot('Position', [0.08+(patchSizeIndex-1)*0.23 0.08 0.18 0.22]);
        csfAxes = subplot('Position', [0.08+(patchSizeIndex-1)*0.23 0.31 0.18 0.68]);
    
        theLegends = {};
        
        hold on;
        
        % Plot the corresponding Crowell&Banks data
        if (~isempty(find(data.patchSize2SigmaCycles(patchSizeIndex)==3.3)))
            plot(csfAxes, CrowellBanksEx3MSB_33cycles.x, CrowellBanksEx3MSB_33cycles.y, 'kv-', 'Color', [0.4 0.4 0.4], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
            'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('MSB, 3.3 cycles');
            
            plot(csfAxes, CrowellBanksEx3JAC_33cycles.x, CrowellBanksEx3JAC_33cycles.y, 'k^-', 'Color', [0.4 0.4 0.4], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
             'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('JAC, 3.3 cycles');
        end
        
        if (~isempty(find(data.patchSize2SigmaCycles(patchSizeIndex)==1.7)))
            plot(csfAxes, CrowellBanksEx3MSB_17cycles.x, CrowellBanksEx3MSB_17cycles.y, 'kv-', 'Color', [0.4 0.4 0.4], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
            'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('MSB, 1.7 cycles');
            
            plot(csfAxes, CrowellBanksEx3JAC_17cycles.x, CrowellBanksEx3JAC_17cycles.y, 'k^-', 'Color', [0.4 0.4 0.4], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
             'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('JAC, 1.7 cycles');
        end
        
        for observerTypeIndex = 1:numel(data.observerTypesExamined)
                
            % ideal observer data
            cpd = squeeze(referenceData.cpd(patchSizeIndex, observerTypeIndex, :));
            contrastSensitivityReference = squeeze(referenceData.contrastSensitivity(patchSizeIndex, observerTypeIndex, :));
            plot(cpd, contrastSensitivityReference, 'o:', 'Color', [1 0 0], ... 
                    'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', markerSize, ...
                    'LineWidth', 1.0); 
            theLegends{numel(theLegends)+1} = referenceData.theLegends{patchSizeIndex};
            
            % computational observer data
            cpd = squeeze(data.cpd(patchSizeIndex, observerTypeIndex, :));
            contrastSensitivity = squeeze(data.contrastSensitivity(patchSizeIndex, observerTypeIndex, :));
            pp = plot(csfAxes, cpd, contrastSensitivity, 'o-', 'Color', [1 0 0], ... 
                    'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', markerSize, ...
                    'LineWidth', 1.5);
   
                
            theLegends{numel(theLegends)+1} = data.theLegends{patchSizeIndex};
            
            plot(ratioAxes, cpd, contrastSensitivity./contrastSensitivityReference, 'o-', 'Color', [1 0 0], ... 
                    'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', markerSize, ...
                    'LineWidth', 1.5);
        end % observerType
        
        for k = 1:numel(theLegends)
            theLegends{k} = strrep(theLegends{k}, 'cycles', 'c');
        end
        theLegend = legend(theLegends);
        set(theLegend, 'FontSize', legendFontSize);
        
        set(csfAxes, 'YScale', 'log', 'XScale', 'log', 'YTick', csTicks, 'XTick', sfTicks, 'XTickLabel', {}, ...
                    'XLim', [sfLims(1)-dx1 sfLims(2)+dx2], 'YLim', [csLims(1)-dy1 csLims(2)+dy2]);
        box(csfAxes, 'on'); grid(csfAxes, 'on');
        set(csfAxes, 'FontSize', fontSize, 'TickLength',[0.02, 0.02], 'LineWidth', axesLineWidth);
        
        if (patchSizeIndex == 1)
            ylabel(csfAxes,'\it contrast sensitivity', 'FontWeight', 'normal', 'FontSize', labelFontSize);
        else
            set(csfAxes, 'YTickLabel', {});
        end
        
        set(csfAxes,'GridLineStyle','-');
        set(csfAxes,'GridAlpha', 0.15);
        set(csfAxes,'MinorGridLineStyle','-')
        set(csfAxes,'MinorGridAlpha', 0.1);

        % The ratios
        set(ratioAxes, 'YScale', 'log', 'XScale', 'log', 'YTick', theRatioTicks, 'YLim', theRatioLims, ...
                          'XTick', sfTicks, 'XLim', [sfLims(1)-dx1 sfLims(2)+dx2]);

        
        
        ytickformat(ratioAxes,'%.2f');
        box(ratioAxes, 'on'); grid(ratioAxes, 'on');
        set(ratioAxes, 'FontSize', fontSize, 'TickLength',[0.02, 0.02] , 'LineWidth', axesLineWidth);

        xlabel(ratioAxes,'\it spatial frequency (c/deg)', 'FontWeight', 'normal', 'FontSize', labelFontSize);
        if (patchSizeIndex == 1)
            ylabel(ratioAxes,'\it sensitivity ratio', 'FontWeight', 'normal', 'FontSize', labelFontSize);
        else
            set(ratioAxes, 'YTickLabel', {});
        end

        set(ratioAxes,'GridLineStyle','-');
        set(ratioAxes,'GridAlpha', 0.15);
        set(ratioAxes,'MinorGridLineStyle','-')
        set(ratioAxes,'MinorGridAlpha', 0.1);

        drawnow;
    
    end % patchSizeIndex

    
    NicePlot.exportFigToPDF(sprintf('CrowellBanks_%s.pdf', performanceClassifier), hFig, 300);
end

