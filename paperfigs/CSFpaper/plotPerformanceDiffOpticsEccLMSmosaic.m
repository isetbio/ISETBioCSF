function plotPerformanceDiffOpticsEccLMSmosaic

    load('/Users/nicolas/Desktop/figDataGeislerOptics.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdGeislerOptics = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesGeislerOptics = 1./(thresholdContrasts*referenceContrast);
    
    load('/Users/nicolas/Desktop/figDataBestPSFSubjectOptics.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdBestPSFSubject= figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesBestPSFSubject = 1./(thresholdContrasts*referenceContrast);
    
    load('/Users/nicolas/Desktop/figDataDefaultSubjectOptics.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdDefaultSubject= figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesDefaultSubject = 1./(thresholdContrasts*referenceContrast);
    
    load('/Users/nicolas/Desktop/figDataTypicalSubjectOptics.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdTypicalSubject= figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesTypicalSubject = 1./(thresholdContrasts*referenceContrast);
    
    % Add the defosued subject
    load('/Users/nicolas/Desktop/figDataDefocusedSubjectOptics.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdDefocusedSubject= figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesDefocusedSubject = 1./(thresholdContrasts*referenceContrast);
    
    % Add the defosued subject
    load('/Users/nicolas/Desktop/figDataVeryDefocusedSubjectOptics.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdVeryDefocusedSubject= figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesVeryDefocusedSubject = 1./(thresholdContrasts*referenceContrast);
    
    
    
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 400 560]);
    
    colors(1,:) = [0.7 0.7 0.7];
    colors(2:6,:) = brewermap(5, 'Set1');
    
    faceColor = [0.1 0.1 0.1];
    
    plot(cpdGeislerOptics, contrastSensitivitiesGeislerOptics, 'o-', 'Color', squeeze(colors(1,:))*0.8, ...
        'MarkerEdgeColor', squeeze(colors(1,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(1,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
    hold on;
    plot(cpdBestPSFSubject, contrastSensitivitiesBestPSFSubject, 'o-', 'Color', squeeze(colors(2,:))*0.8,...
        'MarkerEdgeColor', squeeze(colors(2,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(2,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
    plot(cpdDefaultSubject, contrastSensitivitiesDefaultSubject, 'o-', 'Color', squeeze(colors(3,:))*0.8,...
        'MarkerEdgeColor', squeeze(colors(3,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(3,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
    plot(cpdTypicalSubject, contrastSensitivitiesTypicalSubject, 'o-', 'Color', squeeze(colors(4,:))*0.8, ...
        'MarkerEdgeColor', squeeze(colors(4,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(4,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
   
    plot(cpdDefocusedSubject, contrastSensitivitiesDefocusedSubject, 'o-', 'Color', squeeze(colors(5,:))*0.8,...
        'MarkerEdgeColor', squeeze(colors(5,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(5,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
     plot(cpdVeryDefocusedSubject, contrastSensitivitiesVeryDefocusedSubject, 'o-', 'Color', squeeze(colors(6,:))*0.8,...
        'MarkerEdgeColor', squeeze(colors(6,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(6,:)) + faceColor), ...
        'MarkerSize',12,'LineWidth',2);
    
    
    
    hL = legend({'Geisler PSF', 'best PSF', 'default PSF', 'typical PSF', 'defocused PSF', 'very defocused PSF'});
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
    title('ecc-based LMS mosaic')
    NicePlot.exportFigToPDF('OpticsVsOpticsEccLMSMosaic.pdf', hFig, 300);
    
end


