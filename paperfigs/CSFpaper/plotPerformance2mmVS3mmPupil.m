function plotPerformance2mmVS3mmPupil
    load('/Users/nicolas/Desktop/figData3mmPupil.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpd3mm = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts3mm = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivities3mm = 1./(thresholdContrasts3mm*referenceContrast);
     
    
    load('/Users/nicolas/Desktop/figData2mmPupil.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpd2mm = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts2mm = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivities2mm = 1./(thresholdContrasts2mm*referenceContrast);
    
    hiResSFRange(1) = max([figData.E(:,1); figData.D(:,1); figData.C(:,1)]);
    hiResSFRange(2) = min([figData.E(:,1); figData.D(:,1); figData.C(:,1)]);
    hiResSF = logspace(log10(hiResSFRange(1)), log10(hiResSFRange(2)), 100);
    midLuminanceBanksCurve = interp1(figData.D(:,1), figData.D(:,2), hiResSF);
    
    clear figData
    
    colors= [0 0.7 0];
    
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 400 560]);
    
    plot(hiResSF,midLuminanceBanksCurve,'-','Color', colors,'LineWidth',3);
    hold on;
    plot(cpd2mm, contrastSensitivities2mm, 'ko-', ...
        'MarkerEdgeColor', colors/1.5, ...
        'MarkerFaceColor', min(1, colors + [0.6 0.6 0.6]), ...
        'MarkerSize',11,'LineWidth',2);
    
    colors2 = [1 0.5 0.5];
    plot(cpd3mm, contrastSensitivities3mm, 'ks-', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', colors2, ...
        'MarkerSize',12,'LineWidth',2);
    
    
    hL = legend({'Banks (2mm)', 'cs(2mm)', 'cs(3mm)'});
    set(hL, 'FontName', 'Monaco');
    
    csTicks = [2 5 10 20 50 100 200 500 1000 2000 5000 10000];
    csLims = [2 12000];
    sfTicks = [1 2 5 10 20 50 100];
    sfLims  = [1.2 100];
    dx1 = 0.2;
    dx2 = 20;
    dy1 = 0.4;
    dy2 = 1000;
    
    set(gca, 'YScale', 'log', 'XScale', 'log', 'YTick', csTicks, 'XTick', sfTicks, ...
        'XLim', [sfLims(1)-dx1 sfLims(2)+dx2], 'YLim', [csLims(1)-dy1 csLims(2)+dy2]);
    box on; grid on;
    set(gca, 'FontSize', 14, 'TickLength',[0.02, 0.02]);
    xlabel('spatial frequency (c/deg)', 'FontWeight', 'bold');
    ylabel('contrast sensitivity', 'FontWeight', 'bold');
    NicePlot.exportFigToPDF('Pupil2vs3mmCS.pdf', hFig, 300);
    
    
    hFig = figure(2);
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 400 400]);
    clf;
    hold on;
    
    pupilRatio = 3/2;
    yLim =   [pupilRatio/2 pupilRatio*2];
    yTicks = 0.5:0.5:3;
    
    plot(cpd3mm, contrastSensitivities3mm./contrastSensitivities2mm, 'ko', ...
        'MarkerFaceColor', [0.75 0.75 0.75], ...
        'MarkerSize',11,'LineWidth',2 ...
        );
    plot([1 200], pupilRatio*[1 1], 'k-', 'LineWidth', 1.5);
    set(gca, 'YScale', 'log', 'XScale', 'log', 'YTick', yTicks', 'YTickLabels', sprintf('%2.2f\n',yTicks), 'XTick', sfTicks, ...
        'XLim', [sfLims(1)-dx1 sfLims(2)+dx2], 'YLim', yLim);
    box on; grid on;
    hL = legend({'cs(3mm) : cs(2mm)', '   3mm  : 2mm'});
    set(hL, 'FontName', 'Monaco');
    
    set(gca, 'FontSize', 14, 'TickLength',[0.02, 0.02]);
    xlabel('spatial frequency (c/deg)', 'FontWeight', 'bold');
    ylabel('cs ratio', 'FontWeight', 'bold');
  
    NicePlot.exportFigToPDF('Pupil2vs3mmCSratio.pdf', hFig, 300);
end
