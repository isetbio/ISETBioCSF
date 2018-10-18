function hFig = generateFigureForPaper(theFigData, variedParamLegends, variedParamName, fixedParamName, varargin)

    p = inputParser;
    p.addParameter('figureType', 'CSF', @ischar);
    p.addParameter('figDataIndicesToDisplay', [], @isnumeric);
    p.addParameter('showBanksPaperIOAcurves', false, @islogical);
    p.addParameter('showOnly23CDM2IOAcurve', false, @islogical);
    p.addParameter('showSubjectData', false, @islogical);
    p.addParameter('showSubjectMeanData', false, @islogical);
    p.addParameter('showLegend', true, @islogical);
    p.addParameter('plotFirstConditionInGray', true, @islogical);
    p.addParameter('plotUsingLargeBlueDisks', false, @islogical);
    p.addParameter('plotRatiosOfOtherConditionsToFirst', false, @islogical);
    p.addParameter('theRatioLims', []);
    p.addParameter('theRatioTicks', []);
    p.addParameter('inGraphText', '', @ischar);
    p.addParameter('inGraphTextPos', [], @isnumeric);
    p.addParameter('inGraphTextFontSize', [], @isnumeric);
    p.parse(varargin{:});
    
    figureType = p.Results.figureType;
    figDataIndicesToDisplay = p.Results.figDataIndicesToDisplay;
    if (isempty(figDataIndicesToDisplay))
        figDataIndicesToDisplay = 1:numel(theFigData);
    end
    showBanksPaperIOAcurves = p.Results.showBanksPaperIOAcurves;
    showOnly23CDM2IOAcurve = p.Results.showOnly23CDM2IOAcurve;
    showLegend = p.Results.showLegend;
    showSubjectData = p.Results.showSubjectData;
    showSubjectMeanData = p.Results.showSubjectMeanData;
    plotFirstConditionInGray = p.Results.plotFirstConditionInGray;
    plotUsingLargeBlueDisks = p.Results.plotUsingLargeBlueDisks;
    inGraphText = p.Results.inGraphText;
    inGraphTextPos = p.Results.inGraphTextPos;
    inGraphTextFontSize = p.Results.inGraphTextFontSize;
    plotRatiosOfOtherConditionsToFirst = p.Results.plotRatiosOfOtherConditionsToFirst;
    theRatioLims = p.Results.theRatioLims;
    theRatioTicks = p.Results.theRatioTicks;
    
    if (~isempty(theFigData))
        for condIndex = 1:numel(theFigData)
            figData = theFigData{condIndex};
            lumIndex = 1;
            referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
            cpd{condIndex} = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
            thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
            contrastSensitivity{condIndex} = 1./(thresholdContrasts*referenceContrast);
        end % condIndex
    end
    
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
        load('BanksCSF.mat', 'BanksCSF');
        figData.BanksCSF = BanksCSF;
        
        %banksFactor = 1;
        %plot(theAxes,figData.D(:,1),figData.D(:,2)*banksFactor,'-','Color', squeeze(colors(1,:)), 'LineWidth',2);
        %plot(theAxes,figData.C(:,1),figData.C(:,2)*banksFactor,'-','Color', squeeze(colors(2,:)), 'LineWidth',2);
        %plot(theAxes,figData.E(:,1),figData.E(:,2)*banksFactor,'-','Color', squeeze(colors(3,:)), 'LineWidth',2);

        if (showOnly23CDM2IOAcurve)
            plot(theAxes,figData.BanksCSF('34 cd/m2').x,figData.BanksCSF('34 cd/m2').y,'--','Color', [0 0 0], 'LineWidth',3);
            nLegends = numel(variedParamLegends);
            if (nLegends > 0)
                for kkk = nLegends:-1:1
                  variedParamLegends{kkk+1} = variedParamLegends{kkk};
                end
            end
            variedParamLegends{1} = 'Ideal Observer (Banks et al ''87)';
        else
            plot(theAxes,figData.BanksCSF('34 cd/m2').x,figData.BanksCSF('34 cd/m2').y,'-','Color', squeeze(colors(1,:)), 'LineWidth',2);
            plot(theAxes,figData.BanksCSF('340 cd/m2').x,figData.BanksCSF('340 cd/m2').y,'-','Color', squeeze(colors(2,:)), 'LineWidth',2);
            plot(theAxes,figData.BanksCSF('3.4 cd/m2').x,figData.BanksCSF('3.4 cd/m2').y,'-','Color', squeeze(colors(3,:)), 'LineWidth',2);
        end
        
    end
   
    for condIndex = 1:numel(theFigData)
        
        plotThisCSFAtFullContrast = any(ismember(figDataIndicesToDisplay, condIndex));
        
        if (plotThisCSFAtFullContrast)
            saturationFactor = 0;
        else
            saturationFactor = 0.3;
        end
        
        if (plotUsingLargeBlueDisks)
            if (saturationFactor == 0)
                color1 = [0.3 0.3 1.0];
                color2 = [0 0 1];
                color3 = [0.5 0.5 1.0];
            else
                color1 = [0.9 0.9 1.0];
                color2 = [0.8 0.8 1];
                color3 = [0.9 0.9 1.0];
            end
            plot(theAxes, cpd{condIndex}, contrastSensitivity{condIndex}, 'o-', 'Color', color1, ...
                'MarkerEdgeColor', color2, ...
                'MarkerFaceColor', color3, ...
                'MarkerSize',16,'LineWidth',3);
        else
            if (saturationFactor == 0)
                edgeColor = max(0, squeeze(colors(condIndex,:))-faceColor);
                faceColor2 = min(1, squeeze(colors(condIndex,:)) + faceColor);
                markerSize = 14;
            else
                edgeColor = max(0, squeeze(colors(condIndex,:))-faceColor);
                faceColor2 = min(1, squeeze(colors(condIndex,:)) + faceColor);
                edgeColor  = [0.7 0.7 0.7] + saturationFactor * edgeColor;
                faceColor2 = [0.7 0.7 0.7] + saturationFactor * faceColor2;
                markerSize = 12;
            end
            plot(theAxes, cpd{condIndex}, contrastSensitivity{condIndex}, 'o-', 'Color', edgeColor, ...
                'MarkerEdgeColor', edgeColor, ...
                'MarkerFaceColor', faceColor2, ...
                'MarkerSize',markerSize,'LineWidth',3);
        end
    end
   
    
    if (showSubjectData)
        load('BanksCSF.mat', 'BanksCSF');
        figData.BanksCSF = BanksCSF;
        
        idealObserverData.SFs = figData.BanksCSF('34 cd/m2').x;
        idealObserverData.CSF = figData.BanksCSF('34 cd/m2').y;
        
        msbSubjectSFs =  [ 4.88  6.65  9.8  13.7 19.7  28.0 40];
        msbSubjectCSFs = [72.5  38.4   18.1 14.1  6.82 3.25 1.85];
        
        pjbSubjectSFs =  [4.81 6.81 9.90 13.9 20.0 27.4];
        pjbSubjectCSFs = [61.6 32.6 12.4 8.47 4.96 2.44];
        
        % Compute the mean of the subjectCSF
        [pjbSubjectSFsHiRes, pjbSubjectCSFsHiRes] = fitData(pjbSubjectSFs, pjbSubjectCSFs);
        [msbSubjectSFsHiRes, msbSubjectCSFsHiRes] = fitData(msbSubjectSFs, msbSubjectCSFs);
        
        %plot(pjbSubjectSFsHiRes, pjbSubjectCSFsHiRes, 'r-', 'LineWidth', 1.5);
        %plot(msbSubjectSFsHiRes, msbSubjectCSFsHiRes, 'b-', 'LineWidth', 1.5);
    
        meanSubjectCSF = 0.5*(msbSubjectCSFsHiRes+pjbSubjectCSFsHiRes);
    
        if (~isempty(theFigData))
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
        end
        
        % Add in the subject data
        subjectColor = [0.5 0.5 0.5];
        plot(theAxes, msbSubjectSFs, msbSubjectCSFs, '^', ...
            'MarkerEdgeColor', subjectColor-0.3, ...
            'MarkerFaceColor', subjectColor, ...
            'MarkerSize',14,'LineWidth',2);
        
        subjectColor = [0.8 0.8 0.8];
        plot(theAxes, pjbSubjectSFs, pjbSubjectCSFs, '^', ...
            'MarkerEdgeColor', subjectColor-0.3, ...
            'MarkerFaceColor', subjectColor, ...
            'MarkerSize',14,'LineWidth',2);
        
        
        % Legends
        variedParamLegends{numel(variedParamLegends)+1} = 'Subject MSB (Banks et al ''87)';
        variedParamLegends{numel(variedParamLegends)+1} = 'Subject PJB (Banks et al ''87)';
        
        if (showSubjectMeanData)
            plot(theAxes, msbSubjectSFsHiRes, meanSubjectCSF, 'k-', 'LineWidth', 2);
            variedParamLegends{numel(variedParamLegends)+1} = 'Mean of subjects PJB,MSB';
        end
    end % showSubjectData
    
    % Add legend
    if (showLegend)
        hL = legend(theAxes, variedParamLegends);
    else
        hL = [];
    end
    
    % Add text
    if (~isempty(inGraphText))
        if (isempty(inGraphTextPos))
            if (plotRatiosOfOtherConditionsToFirst)
                inGraphTextPos(1) = 1.6;
                inGraphTextPos(2) = 7000;
            else
                inGraphTextPos(1) = 1.6;
                inGraphTextPos(2) = 7000;
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
            
            saturationFactor = 0;
            edgeColor = max(0, squeeze(colors(condIndex,:))-faceColor);
            faceColor2 = min(1, squeeze(colors(condIndex,:)) + faceColor);
                
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
                'ko-', 'Color', edgeColor, ...
                'MarkerEdgeColor', edgeColor, ...
                'MarkerFaceColor', faceColor2, ...
                'MarkerSize',14,'LineWidth',3);
        end
        
        if (showSubjectData)
            plot(theRatioAxes, refCPD, meanSubjectRatios, ...
                'k-', 'Color', [0 0 0], ...
                'MarkerEdgeColor', [0 0 0], ...
                'MarkerFaceColor', [0.5 0.5 0.5], ...
                'MarkerSize',14,'LineWidth',3);
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