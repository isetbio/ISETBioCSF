function plotPerformanceDiffClassifiers

    load('/Users/nicolas/Desktop/figDataNOEMMLPT2.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdMLPT = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesMLPT = 1./(thresholdContrasts*referenceContrast);
    
    load('/Users/nicolas/Desktop/figDataBestPSFSubjectOpticsSVM.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdSVM= figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesSVM = 1./(thresholdContrasts*referenceContrast);
    
    load('/Users/nicolas/Desktop/figDataBestPSFSubjectOpticsSVMV1Filter.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdSVMV1filter= figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesSVMV1filter = 1./(thresholdContrasts*referenceContrast);
    
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 400 560]);
    
    faceColor = [0.1 0.1 0.1];
    colors = [0.5 0.5 0.5; 0.8 0.0 0.0; .2 0.7 0.1 ];
    plot(cpdMLPT, contrastSensitivitiesMLPT, 'o-', 'Color', squeeze(colors(1,:))*0.8, ...
        'MarkerEdgeColor', squeeze(colors(1,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(1,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
    hold on;
    plot(cpdSVM, contrastSensitivitiesSVM, 'o-', 'Color', squeeze(colors(2,:))*0.8, ...
        'MarkerEdgeColor', squeeze(colors(2,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(2,:)) + 2*faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
    plot(cpdSVMV1filter, contrastSensitivitiesSVMV1filter, 'o-', 'Color', squeeze(colors(3,:))*0.8, ...
        'MarkerEdgeColor', squeeze(colors(3,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(3,:)) + 2*faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
    hL = legend({'MLPT', 'SVM', 'SVM V1filter'});
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
    title('best PSF subject optics')
    NicePlot.exportFigToPDF('ClassiferPerformance.pdf', hFig, 300);
    
end


