function plotPerformanceEyeMovements

    load('/Users/nicolas/Desktop/figDataNOEMMLPT2.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdNOEMMLPT = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesNOEMMLPT = 1./(thresholdContrasts*referenceContrast);
    
    load('/Users/nicolas/Desktop/figDataNOEMSVMV1Filter2.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdNOEMSVMV1Filter = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesNOEMSVMV1Filter = 1./(thresholdContrasts*referenceContrast);
    
    load('/Users/nicolas/Desktop/figDataNOEMSVM.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdNOEMSVM = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesNOEMSVM = 1./(thresholdContrasts*referenceContrast);
    
    load('/Users/nicolas/Desktop/figDataEMSVM.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdEMSVM = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesEMSVM = 1./(thresholdContrasts*referenceContrast);
    
    load('/Users/nicolas/Desktop/figDataEMSVMV1filter.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdEMSVMV1filter = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesEMSVMV1filter = 1./(thresholdContrasts*referenceContrast);
    
    
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 400 560]);
    
    colors(1,:) = [0.5 0.5 0.5];
    colors(2:5,:) = brewermap(4, 'Paired');
    faceColor = [0.1 0.1 0.1];
    
    plot(cpdNOEMMLPT, contrastSensitivitiesNOEMMLPT, 'o-', 'Color', squeeze(colors(1,:))*0.8, ...
        'MarkerEdgeColor', squeeze(colors(1,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(1,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
    hold on;
    plot(cpdNOEMSVM, contrastSensitivitiesNOEMSVM, 'o-', 'Color', squeeze(colors(2,:))*0.8, ...
        'MarkerEdgeColor', squeeze(colors(2,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(2,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
    plot(cpdNOEMSVMV1Filter, contrastSensitivitiesNOEMSVMV1Filter, 'o-', 'Color', squeeze(colors(3,:))*0.8, ...
        'MarkerEdgeColor', squeeze(colors(3,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(3,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
    plot(cpdEMSVM , contrastSensitivitiesEMSVM, 'o-', 'Color', squeeze(colors(4,:))*0.8, ...
        'MarkerEdgeColor', squeeze(colors(4,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(4,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
    plot(cpdEMSVMV1filter, contrastSensitivitiesEMSVMV1filter, 'o-', 'Color', squeeze(colors(5,:))*0.8, ...
        'MarkerEdgeColor', squeeze(colors(5,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(5,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
    hL = legend({'no em (MLPT)', 'no em (SVM)', 'no em (SVM-V1)', 'em (SVM)', 'em (SVM-V1)'});
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
    NicePlot.exportFigToPDF('PerformanceEyeMovement.pdf', hFig, 300);
    
end


