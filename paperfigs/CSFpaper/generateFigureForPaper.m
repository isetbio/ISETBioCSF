function generateFigureForPaper(theFigData, variedParamLegends, variedParamName, fixedParamName, varargin)

    p = inputParser;
    p.addParameter('figureType', 'CSF', @ischar);
    p.addParameter('showBanksPaperIOAcurves', false, @islogical);
    p.addParameter('showSubjectData', false, @islogical);
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
    showSubjectData = p.Results.showSubjectData;
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
        if (numel(theFigData)>1) || (showSubjectData)
            colors(2:numel(theFigData),:) = brewermap(numel(theFigData)-1, 'Set1');
            if (showSubjectData)
                colorsNum = numel(theFigData)+2;
            else
                colorsNum = numel(theFigData);
            end
            colors(2:colorsNum+1,:) = brewermap(colorsNum, 'Set1');
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
    
    if (showSubjectData)
        load('BanksCSF.mat', 'BanksCSF');
        figData.BanksCSF = BanksCSF;
        
        idealObserverData.SFs = figData.BanksCSF('34 cd/m2').x;
        idealObserverData.CSF = figData.BanksCSF('34 cd/m2').y;
        msbSubjectRatios = [...
            5.0910, 16.0268
            7.0576,   18.3709
            10.0095,   21.7148
            14.1960,   17.7137
            20.4416,   16.6040
            28.5542,   16.0562
            40.8058,    9.1792];
        
        pjbSubjectRatios = [...
            4.92,	17.57
            6.93,	21.20
            8.89,	31.81
            13.86,	29.48
            19.59,	21.29
            27.87,	18.54];

        msbSubjectSFs = squeeze(msbSubjectRatios(:,1));
        msbSubjectRatios = squeeze(msbSubjectRatios(:,2));
        
        pjbSubjectSFs = squeeze(pjbSubjectRatios(:,1));
        pjbSubjectRatios = squeeze(pjbSubjectRatios(:,2));
        
        for dataPoint = 1:numel(msbSubjectSFs)
            targetSF = msbSubjectSFs(dataPoint);
            [~,idx] = min(abs(targetSF-idealObserverData.SFs));
            %fprintf('matched SFs: target:%2.2f ideal observer match:%2.2f\n', targetSF, idealObserverData.SFs(idx));
            msbSubjectCSFs(dataPoint) = idealObserverData.CSF(idx)./msbSubjectRatios(dataPoint);
        end

        for dataPoint = 1:numel(pjbSubjectRatios)
            targetSF = pjbSubjectSFs(dataPoint);
            [~,idx] = min(abs(targetSF-idealObserverData.SFs));
            fprintf('matched SFs: target:%2.2f ideal observer match:%2.2f\n', targetSF, idealObserverData.SFs(idx));
            pjbSubjectCSFs(dataPoint) = idealObserverData.CSF(idx)./pjbSubjectRatios(dataPoint);
        end
        
        % Compute the mean of the subjectCSF
        [pjbSubjectSFsHiRes, pjbSubjectCSFsHiRes] = fitData(pjbSubjectSFs, pjbSubjectCSFs);
        plot(pjbSubjectSFsHiRes, pjbSubjectCSFsHiRes, 'r-', 'LineWidth', 1.5);
    
        [msbSubjectSFsHiRes, msbSubjectCSFsHiRes] = fitData(msbSubjectSFs, msbSubjectCSFs);
        plot(msbSubjectSFsHiRes, msbSubjectCSFsHiRes, 'b-', 'LineWidth', 1.5);
    
        meanSubjectCSF = 0.5*(msbSubjectCSFsHiRes+pjbSubjectCSFsHiRes);
    
        subjectCPD = msbSubjectSFsHiRes;
        refSensitivity = contrastSensitivity{1};
        refCPD = cpd{1};
        
        for k = 1:numel(refCPD)
            idx = find(msbSubjectSFsHiRes == refCPD(k));
            if (~isempty(idx))
                meanSubjectRatios(k) = meanSubjectCSF(idx)./refSensitivity(k);
            else
                meanSubjectRatios(k) = nan;
            end
        end
        
        % Add in the subject data
        colorIndex = size(colors,1)-1;
        plot(theAxes, msbSubjectSFs, msbSubjectCSFs, '^', 'Color', squeeze(colors(colorIndex,:)), ...
            'MarkerEdgeColor', max(0, squeeze(colors(colorIndex,:))-faceColor), ...
            'MarkerFaceColor', min(1, squeeze(colors(colorIndex,:)) + faceColor), ...
            'MarkerSize',12,'LineWidth',2);
        
        colorIndex = size(colors,1);
        plot(theAxes, pjbSubjectSFs, pjbSubjectCSFs, '^', 'Color', squeeze(colors(colorIndex,:)), ...
            'MarkerEdgeColor', max(0, squeeze(colors(colorIndex,:))-faceColor), ...
            'MarkerFaceColor', min(1, squeeze(colors(colorIndex,:)) + faceColor), ...
            'MarkerSize',12,'LineWidth',2);
        
        plot(theAxes, msbSubjectSFsHiRes, meanSubjectCSF, 'k-', 'LineWidth', 2);
        
        variedParamLegends{numel(variedParamLegends)+1} = 'Subject MSB (Banks et al ''87)';
        variedParamLegends{numel(variedParamLegends)+1} = 'Subject PJB (Banks et al ''87)';
        variedParamLegends{numel(variedParamLegends)+1} = 'Mean of subjects PJB,MSB';
   end
    
    % Add legend
    hL = legend(theAxes, variedParamLegends);
    
    % Add text
    if (~isempty(inGraphText))
        if (isempty(inGraphTextPos))
            if (plotRatiosOfOtherConditionsToFirst)
                inGraphTextPos(1) = 1.6;
                inGraphTextPos(2) = 10000;
            else
                inGraphTextPos(1) = 1.6;
                inGraphTextPos(2) = 9500;
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
        
        if (showSubjectData)
            plot(theRatioAxes, refCPD, meanSubjectRatios, ...
                'k-', 'Color', [0 0 0], ...
                'MarkerEdgeColor', [0 0 0], ...
                'MarkerFaceColor', [0.5 0.5 0.5], ...
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
    
    %if (strcmp(fixedParamName, 'ComparedToBanks87Photocurrents'))
    %    set(hL, 'Location', 'NorthOutside');
    %end
    
    exportsDir = strrep(isetRootPath(), 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    figureName = fullfile(exportsDir, sprintf('%sVary%s.pdf', variedParamName, fixedParamName));
    NicePlot.exportFigToPDF(figureName, hFig, 300);
end

function [xFit, yFit] = fitData(xData,yData)

    xdata = xData(:);
    ydata = yData(:);
    
    F = @(p,xdata)p(1) + p(2)*exp(-p(3)*xdata.^p(6)) + p(4)*exp(-p(5)*xdata.^2);
    params0 = [1.3090   43.5452    0.0379   92.5870    0.0356    1.3392]; 
    
    [params,resnorm,~,exitflag,output] = lsqcurvefit(F,params0,xdata,ydata);
    
    xFit = 4:0.1:32;
    yFit = F(params,xFit);

end