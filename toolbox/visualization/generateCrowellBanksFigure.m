function generateCrowellBanksFigure(data)
    
    performanceClassifier = 'svmV1FilterEnsemble';
    %performanceClassifier = 'svmV1FilterBank';
    
    
    if notDefined('data')
        
        dataFilename = sprintf('data_mlpt5_0_4.mat'); 
        load(fullfile('/Users/nicolas/Desktop',dataFilename), 'data');
        reference_data5_0_4 = data;
        
        dataFilename = sprintf('data_%s5_0_4.mat', strrep(performanceClassifier, ' ', '_'));
        load(fullfile('/Users/nicolas/Desktop',dataFilename), 'data');
        data5_0_4 = data;
        
        
        dataFilename = sprintf('data_mlpt.mat');
        load(fullfile('/Users/nicolas/Desktop',dataFilename), 'data');
        referenceData = data;
        
        
        dataFilename = sprintf('data_%s.mat', strrep(performanceClassifier, ' ', '_'));
        load(fullfile('/Users/nicolas/Desktop',dataFilename), 'data');
        
    end
    
        
    
    csTicks = [2 5 10 20 50 100 200 500 1000 2000 5000 10000];
    csLims = [1.5 10000];
    sfTicks = [1 2 5 10 20 50 100];
    sfLims  = [1.5 80];
    dx1 = 0.2; dx2 = 10; dy1 = 0.5; dy2 = 1000;
    theRatioTicks = [0.01 0.03 0.1 0.3 1];
    theRatioLims = [0.01 0.5];
    
    hFig = figure(1234); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 800 390]);
    
    
    
    markerSize = 6;
    lineWidth = 1.15;
    fontSize = 12;
    legendFontSize = 7;
    labelFontSize = 12;
    axesLineWidth = 0.25;
           
    % Load CrowellBanks data set
    load('CrowellBanksSubjects.mat');
    CrowellBanksIdealObserver_66cycles.x = [...
        4.718459307587E0 ...
        7.260488284366E0 ...
        1.154494509231E1 ...
        1.746534887941E1 ...
        2.513286550973E1 ...
        3.220662661559E1 ...
        3.925254462762E1 ...
        4.704464279244E1 ...
        5.273168658502E1 ...
        ];
    
    CrowellBanksIdealObserver_66cycles.y = [...
        2.919940226722E3 ...
        1.708783259142E3 ...
        8.417977415172E2 ...
        4.227063589600E2 ...
        2.004190772077E2 ...
        1.128837891685E2 ...
        5.668427598655E1 ...
        2.687589160797E1 ...
        1.052351113242E1 ];
        

    % Plot the computed data
    for patchSizeIndex  = 1:4
        
        ratioAxes = subplot('Position', [0.08+(patchSizeIndex-1)*0.21 0.09 0.20 0.22]);
        csfAxes = subplot('Position', [0.08+(patchSizeIndex-1)*0.21 0.32 0.20 0.67]);
    
        theLegends = {};
        
        hold on;
        
        % Plot the Crowell&Banks data
        if (patchSizeIndex ==  1)
            
            % Crowell&Banks ideal observer data
            plot(CrowellBanksIdealObserver_66cycles.x, CrowellBanksIdealObserver_66cycles.y*3.3/6.6, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 2.0);
            theLegends{numel(theLegends)+1} = sprintf('C&B ideal observer, 3.3 cycles');
            
            % Plot the Crowell&Banks psychophysical data for subject MSB
            plot(csfAxes, CrowellBanksEx3MSB_33cycles.x, CrowellBanksEx3MSB_33cycles.y, 'kv-', 'Color', [0.4 0.4 0.4], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
            'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('MSB, 3.3 cycles');
            
            % Plot the Crowell&Banks psychophysical data for subject JAC
            plot(csfAxes, CrowellBanksEx3JAC_33cycles.x, CrowellBanksEx3JAC_33cycles.y, 'k^-', 'Color', [0.4 0.4 0.4], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
             'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('JAC, 3.3 cycles');
        end
        
        if (patchSizeIndex == 2)
            
            % Crowell&Banks ideal observer data
            plot(CrowellBanksIdealObserver_66cycles.x, CrowellBanksIdealObserver_66cycles.y*1.7/6.6, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 2.0);
            theLegends{numel(theLegends)+1} = sprintf('C&B ideal observer, 1.7 cycles');
            
            % Plot the Crowell&Banks psychophysical data for subject MSB
            plot(csfAxes, CrowellBanksEx3MSB_17cycles.x, CrowellBanksEx3MSB_17cycles.y, 'kv-', 'Color', [0.4 0.4 0.4], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
            'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('MSB, 1.7 cycles');
            
            % Plot the Crowell&Banks psychophysical data for subject JAC
            plot(csfAxes, CrowellBanksEx3JAC_17cycles.x, CrowellBanksEx3JAC_17cycles.y, 'k^-', 'Color', [0.4 0.4 0.4], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
             'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('JAC, 1.7 cycles');
        end
        
        if (patchSizeIndex == 3)
            
            % Crowell&Banks ideal observer data
            plot(CrowellBanksIdealObserver_66cycles.x, CrowellBanksIdealObserver_66cycles.y*0.8/6.6, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 2.0);
            theLegends{numel(theLegends)+1} = sprintf('C&B ideal observer, 0.8 cycles');
            
            % Plot the Crowell&Banks psychophysical data for subject MSB
            plot(csfAxes, CrowellBanksEx3MSB_08cycles.x, CrowellBanksEx3MSB_08cycles.y, 'kv-', 'Color', [0.4 0.4 0.4], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
            'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('MSB, 0.8 cycles');
            
            % Plot the Crowell&Banks psychophysical data for subject JAC
            plot(csfAxes, CrowellBanksEx3JAC_08cycles.x, CrowellBanksEx3JAC_08cycles.y, 'k^-', 'Color', [0.4 0.4 0.4], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
             'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('JAC, 0.8 cycles');
        end
        
        if (patchSizeIndex == 4)
            
            % Crowell&Banks ideal observer data
            plot(CrowellBanksIdealObserver_66cycles.x, CrowellBanksIdealObserver_66cycles.y*0.4/6.6, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 2.0);
            theLegends{numel(theLegends)+1} = sprintf('C&B ideal observer, 0.4 cycles');
            
            % Plot the Crowell&Banks psychophysical data for subject MSB
            plot(csfAxes, CrowellBanksEx3MSB_04cycles.x, CrowellBanksEx3MSB_04cycles.y, 'kv-', 'Color', [0.4 0.4 0.4], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
            'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('MSB, 0.4 cycles');
            
            % Plot the Crowell&Banks psychophysical data for subject JAC
            plot(csfAxes, CrowellBanksEx3JAC_04cycles.x, CrowellBanksEx3JAC_04cycles.y, 'k^-', 'Color', [0.4 0.4 0.4], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
             'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('JAC, 0.4 cycles');
        end
        
        
        for observerTypeIndex = 1:numel(data.observerTypesExamined)
                
            % ideal observer data
            if (~isempty(referenceData)) && (size(referenceData.cpd,1) >= patchSizeIndex)
                cpdReference = squeeze(referenceData.cpd(patchSizeIndex, observerTypeIndex, :));
                contrastSensitivityReference = squeeze(referenceData.contrastSensitivity(patchSizeIndex, observerTypeIndex, :));
                plot(cpdReference, contrastSensitivityReference, 'o:', 'Color', [1 0 0], ... 
                    'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', markerSize, ...
                    'LineWidth', 2.0); 
                theLegends{numel(theLegends)+1} = referenceData.theLegends{patchSizeIndex};
            end
            
            if (size(data.cpd,1) >= patchSizeIndex) 
                % computational observer data
                cpd = squeeze(data.cpd(patchSizeIndex, observerTypeIndex, :));
                contrastSensitivity = squeeze(data.contrastSensitivity(patchSizeIndex, observerTypeIndex, :));
                plot(csfAxes, cpd, contrastSensitivity, 'o-', 'Color', [1 0 0], ... 
                        'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', markerSize, ...
                        'LineWidth', 1.5);
            
                theLegends{numel(theLegends)+1} = data.theLegends{patchSizeIndex};
                %if (~isempty(referenceData)) && (size(referenceData.cpd,1) >= patchSizeIndex)
                   ratios = [];
                   cpdRatios = [];
                   for k = 1:numel(cpdReference)
                       idx = find(cpd == cpdReference(k));
                       if (~isempty(idx))
                           cpdRatios(numel(cpdRatios)+1) = cpdReference(k);
                           ratios(numel(ratios)+1) = contrastSensitivity(idx)./contrastSensitivityReference(k);
                       end
                   end
                   plot(ratioAxes, cpdRatios, ratios, 'o-', 'Color', [1 0 0], ... 
                            'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', markerSize, ...
                            'LineWidth', 1.5);
                %end
            end
       
            if (patchSizeIndex == 4)
                cpdReference = squeeze(reference_data5_0_4.cpd(1, observerTypeIndex, :));
                contrastSensitivityReference = squeeze(reference_data5_0_4.contrastSensitivity(1, observerTypeIndex, :));
                plot(cpdReference, contrastSensitivityReference, 'o:', 'Color', [1 0 0], ... 
                    'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', markerSize, ...
                    'LineWidth', 2.0); 
                theLegends{numel(theLegends)+1} = reference_data5_0_4.theLegends{1};
                
                cpd = squeeze(data5_0_4.cpd(1, observerTypeIndex, :));
                contrastSensitivity = squeeze(data5_0_4.contrastSensitivity(1, observerTypeIndex, :));
                plot(csfAxes, cpd, contrastSensitivity, 'o-', 'Color', [1 0 0], ... 
                        'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', markerSize, ...
                        'LineWidth', 1.5);
                theLegends{numel(theLegends)+1} = data5_0_4.theLegends{1};    
            end
            
    
        end % observerType
        
        for k = 1:numel(theLegends)
            theLegends{k} = strrep(theLegends{k}, 'cycles', 'c');
        end
        theLegend = legend(theLegends);
        set(theLegend, 'FontSize', legendFontSize);
        legendPos = get(theLegend, 'Position');
        legendPos(2) = 0.86;
        set(theLegend, 'Position', legendPos);
        
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

