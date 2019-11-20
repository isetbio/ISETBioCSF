function generateCrowellBanksFigureTest(data)
    
    performanceClassifier = 'svmV1FilterEnsemble';
    performanceClassifier2 = 'svmV1FilterBank';
    
    
    
    if notDefined('data')
   
        referenceData = [];
        dataFilename = sprintf('data_mlpt.mat');
        load(fullfile('/Users/nicolas/Desktop',dataFilename), 'data');
        referenceData = data;
        
        data2 = [];
        dataFilename = sprintf('data_%s.mat', strrep(performanceClassifier2, ' ', '_'));
        load(fullfile('/Users/nicolas/Desktop',dataFilename), 'data');
        data2 = data;
        data2.legend = 'energy SVM (single)';
        
        data = [];
        dataFilename = sprintf('data_%s.mat', strrep(performanceClassifier, ' ', '_'));
        load(fullfile('/Users/nicolas/Desktop',dataFilename), 'data');
        data.legend = 'energy SVM (ensemble)';
        
        
    end
    
        
    
    csTicks = [2 5 10 20 50 100 200 500 1000 2000 5000 10000];
    csLims = [1.5 14000];
    sfTicks = [1 2 5 10 20 50 100];
    sfLims  = [1.5 80];
    dx1 = 0.2; dx2 = 10; dy1 = 0.5; dy2 = 1000;
    theRatioTicks = [0.01 0.03 0.1 0.3 1];
    theRatioLims = [0.01 0.5];
    
    hFig = figure(1234); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 800 440]);
    
    
    
    markerSize = 8;
    lineWidth = 1.5;
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
        

    % Plot the data
    availablePatchSizes = [3.3 1.7 0.8 0.4];
    
    for patchSizeIndex = 1:numel(availablePatchSizes)
        
        ratioAxes = subplot('Position', [0.08+(patchSizeIndex-1)*0.21 0.09 0.20 0.22]);
        csfAxes = subplot('Position', [0.08+(patchSizeIndex-1)*0.21 0.32 0.20 0.64]);
    
        theLegends = {};
        
        hold on;
        
        patchSize = availablePatchSizes(patchSizeIndex);
        
        % Plot the Crowell&Banks data
        switch (patchSize)
            case 3.3
                MSB_sensitivity.x = CrowellBanksEx3MSB_33cycles.x;
                MSB_sensitivity.y = CrowellBanksEx3MSB_33cycles.y;
                JAC_sensitivity.x = CrowellBanksEx3JAC_33cycles.x;
                JAC_sensitivity.y = CrowellBanksEx3JAC_33cycles.y;
            case 1.7
                MSB_sensitivity.x = CrowellBanksEx3MSB_17cycles.x;
                MSB_sensitivity.y = CrowellBanksEx3MSB_17cycles.y;
                JAC_sensitivity.x = CrowellBanksEx3JAC_17cycles.x;
                JAC_sensitivity.y = CrowellBanksEx3JAC_17cycles.y;
            case 0.8
                MSB_sensitivity.x = CrowellBanksEx3MSB_08cycles.x;
                MSB_sensitivity.y = CrowellBanksEx3MSB_08cycles.y;
                JAC_sensitivity.x = CrowellBanksEx3JAC_08cycles.x;
                JAC_sensitivity.y = CrowellBanksEx3JAC_08cycles.y;
            case 0.4
                MSB_sensitivity.x = CrowellBanksEx3MSB_04cycles.x;
                MSB_sensitivity.y = CrowellBanksEx3MSB_04cycles.y;
                JAC_sensitivity.x = CrowellBanksEx3JAC_04cycles.x;
                JAC_sensitivity.y = CrowellBanksEx3JAC_04cycles.y;
            otherwise
                error('No sensitivities for patch size %f\n', patchSize);
        end
            
        
        % Crowell&Banks ideal observer data
        plot(CrowellBanksIdealObserver_66cycles.x, CrowellBanksIdealObserver_66cycles.y*availablePatchSizes(patchSizeIndex)/6.6, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2.0);
        theLegends{numel(theLegends)+1} = sprintf('C&B ideal observer');
            
        % Plot the Crowell&Banks psychophysical data for subject MSB
        plot(csfAxes, MSB_sensitivity.x , MSB_sensitivity.y, 'kv-', 'Color', [0.3 0.3 0.5], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
            'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
            theLegends{numel(theLegends)+1} = sprintf('subject MSB');
            
        % Plot the Crowell&Banks psychophysical data for subject JAC
        plot(csfAxes, JAC_sensitivity.x, JAC_sensitivity.y, 'k^-', 'Color', [0.5 0.3 0.3], ...
            'MarkerFaceColor', 0.7*[1 1 1], ...
             'MarkerSize', markerSize*1.4, 'LineWidth', lineWidth);
        theLegends{numel(theLegends)+1} = sprintf('subject JAC');
        
        
        observerTypeIndex = 1;
       
        % plot the ideal observer data
        if (~isempty(referenceData)) && (~isempty(find(referenceData.patchSize2SigmaCycles == patchSize)))
            theCorrespondingPatchSizeIndex = find(referenceData.patchSize2SigmaCycles == patchSize);
            cpdReference = squeeze(referenceData.cpd(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
            contrastSensitivityReference = squeeze(referenceData.contrastSensitivity(theCorrespondingPatchSizeIndex, observerTypeIndex, :));

            validIndices = find(~isnan(cpdReference));
            plot(cpdReference(validIndices ), contrastSensitivityReference(validIndices ), ':', 'Color', [1 0 0], ... 
                    'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', markerSize, ...
                    'LineWidth', 2.0); 
            theLegends{numel(theLegends)+1} = sprintf('ideal observer');
        end
        
        %  plot the compuitational observer data
        if (~isempty(data)) && (~isempty(find(data.patchSize2SigmaCycles == patchSize)))
            theCorrespondingPatchSizeIndex = find(data.patchSize2SigmaCycles == patchSize);
            cpd = squeeze(data.cpd(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
            contrastSensitivity= squeeze(data.contrastSensitivity(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
            plot(cpd, contrastSensitivity, 'o-', 'Color', [1 0 0], ... 
                    'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', markerSize, ...
                    'LineWidth', 2.0); 
            theLegends{numel(theLegends)+1} = sprintf('%s', data.legend);
        end

        %  plot the compuitational observer data2
        if (~isempty(data2)) && (~isempty(find(data2.patchSize2SigmaCycles == patchSize)))
            theCorrespondingPatchSizeIndex = find(data2.patchSize2SigmaCycles == patchSize);
            cpd = squeeze(data2.cpd(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
            contrastSensitivity= squeeze(data2.contrastSensitivity(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
            plot(cpd, contrastSensitivity, 'o-', 'Color', [0 0 1], ... 
                    'MarkerFaceColor', [0.5 0.5 1], 'MarkerSize', markerSize, ...
                    'LineWidth', 2.0); 
            theLegends{numel(theLegends)+1} = sprintf('%s', data2.legend);
        end
        title(sprintf('%2.1f cycles', availablePatchSizes(patchSizeIndex)));
        
        % Plot the ratios
        if (~isempty(referenceData)) && (~isempty(data)) && ...
           (~isempty(find(referenceData.patchSize2SigmaCycles == patchSize))) && ...
           (~isempty(find(data.patchSize2SigmaCycles == patchSize))) 
              
            theCorrespondingPatchSizeIndex = find(referenceData.patchSize2SigmaCycles == patchSize);
            cpdReference = squeeze(referenceData.cpd(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
            contrastSensitivityReference = squeeze(referenceData.contrastSensitivity(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
            
            theCorrespondingPatchSizeIndex = find(data.patchSize2SigmaCycles == patchSize);
            cpd = squeeze(data.cpd(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
            contrastSensitivity = squeeze(data.contrastSensitivity(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
 
            plot(ratioAxes, cpd, contrastSensitivity./contrastSensitivityReference, 'o-', 'Color', [1 0 0], ... 
                    'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', markerSize, ...
                    'LineWidth', 2.0); 
        end
        
        if (~isempty(referenceData)) && (~isempty(data2)) && ...
           (~isempty(find(referenceData.patchSize2SigmaCycles == patchSize))) && ...
           (~isempty(find(data2.patchSize2SigmaCycles == patchSize))) 
              
            theCorrespondingPatchSizeIndex = find(referenceData.patchSize2SigmaCycles == patchSize);
            cpdReference = squeeze(referenceData.cpd(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
            contrastSensitivityReference = squeeze(referenceData.contrastSensitivity(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
            
            theCorrespondingPatchSizeIndex = find(data2.patchSize2SigmaCycles == patchSize);
            cpd = squeeze(data2.cpd(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
            contrastSensitivity = squeeze(data2.contrastSensitivity(theCorrespondingPatchSizeIndex, observerTypeIndex, :));
            hold(ratioAxes, 'on');
            plot(ratioAxes, cpd, contrastSensitivity./contrastSensitivityReference, 'o-', 'Color', [0 0 1 ], ... 
                    'MarkerFaceColor', [0.5 0.5 1], 'MarkerSize', markerSize, ...
                    'LineWidth', 2.0); 
            
        end
        
        for k = 1:numel(theLegends)
            theLegends{k} = strrep(theLegends{k}, 'cycles', 'c');
        end
        
        theLegend = legend(theLegends);
        set(theLegend, 'FontSize', legendFontSize);
        legendPos = get(theLegend, 'Position');
        legendPos(1) = legendPos(1)+0.004;
        legendPos(2) = 0.825;
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
