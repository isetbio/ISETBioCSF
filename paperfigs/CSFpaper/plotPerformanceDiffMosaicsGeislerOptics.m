function plotPerformanceDiffMosaicsGeislerOptics

    load('/Users/nicolas/Desktop/figDataGeislerOpticsBanksMosaic.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdBanksMosaic = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesBanksMosaic = 1./(thresholdContrasts*referenceContrast);
    
    load('/Users/nicolas/Desktop/figDataGeislerOpticsEccLMSMosaic.mat', 'figData');
    lumIndex = 1;
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    cpdEccLMSMosaic= figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
    thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
    contrastSensitivitiesEccLMSMosaic = 1./(thresholdContrasts*referenceContrast);
    
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 400 560]);
    
    colors= [0.0 0.5 1.0; 1 0.0 0.3];
    plot(cpdBanksMosaic, contrastSensitivitiesBanksMosaic, 'bo-', ...
        'MarkerEdgeColor', squeeze(colors(1,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(1,:)) + [0.5 0.5 0.5]), ...
        'MarkerSize',12,'LineWidth',2);
    
    hold on;
    plot(cpdEccLMSMosaic, contrastSensitivitiesEccLMSMosaic, 'ro-', ...
        'MarkerEdgeColor', squeeze(colors(2,:)), ...
        'MarkerFaceColor', min(1, squeeze(colors(2,:)) + [0.5 0.5 0.5]), ...
        'MarkerSize',12,'LineWidth',2);
    
    hL = legend({'Banks LM mosaic', 'ecc-based LMS mosaic'});
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
    title('Geisler optics')
    NicePlot.exportFigToPDF('BanksMosaicVsEccLMSMosaic.pdf', hFig, 300);
    
end


