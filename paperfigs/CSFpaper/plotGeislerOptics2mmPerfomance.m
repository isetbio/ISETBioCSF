function plotGeislerOptics2mmPerfomance
    load('/Users/nicolas/Desktop/figData3mmPupil.mat', 'figData')
    
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 400 560]);
    
    luminances = [3.4 34 340];
    nLuminances = numel(luminances);
    
    referenceContrast = figData.banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1);
    
    colors = [1 0 0; 0 0.7 0; 0 0 1];
    
    csTicks = [2 5 10 20 50 100 200 500 1000 2000 5000 10000];
    csLims = [2 12000];
    sfTicks = [1 2 5 10 20 50 100];
    sfLims  = [1.2 100];
    dx1 = 0.2;
    dx2 = 20;
    dy1 = 0.4;
    dy2 = 1000;
 
    
    % Unshifted banks data
    hiResSFRange(1) = max([figData.E(:,1); figData.D(:,1); figData.C(:,1)]);
    hiResSFRange(2) = min([figData.E(:,1); figData.D(:,1); figData.C(:,1)]);
    hiResSF = logspace(log10(hiResSFRange(1)), log10(hiResSFRange(2)), 100);
    
    minLuminanceBanksCurve = interp1(figData.E(:,1), figData.E(:,2), hiResSF);
    midLuminanceBanksCurve = interp1(figData.D(:,1), figData.D(:,2), hiResSF);
    maxLuminanceBanksCurve = interp1(figData.C(:,1), figData.C(:,2), hiResSF);
    
    
    plot(hiResSF,minLuminanceBanksCurve,'-','Color', colors(1,:),'LineWidth',2);
    hold on
    plot(hiResSF,midLuminanceBanksCurve,'-','Color', colors(2,:),'LineWidth',2);
    plot(hiResSF,maxLuminanceBanksCurve,'-','Color', colors(3,:),'LineWidth',2);
    
    plot(hiResSF,minLuminanceBanksCurve,'k-','Color', min(1, squeeze(colors(1,:)+ [0.6 0.6 0.6])),'LineWidth',3);
    plot(hiResSF,midLuminanceBanksCurve,'k-','Color', min(1, squeeze(colors(2,:)+ [0.6 0.6 0.6])),'LineWidth',4);
    plot(hiResSF,maxLuminanceBanksCurve,'k-','Color', min(1, squeeze(colors(3,:)+ [0.6 0.6 0.6])),'LineWidth',3);
    plot(hiResSF,minLuminanceBanksCurve,'k-','Color', squeeze(colors(1,:)),'LineWidth',1.5);
    plot(hiResSF,midLuminanceBanksCurve,'k-','Color', squeeze(colors(2,:)),'LineWidth',2);
    plot(hiResSF,maxLuminanceBanksCurve,'k-','Color', squeeze(colors(3,:)),'LineWidth',1.5);
    
    banksSquareRootLaw1 = midLuminanceBanksCurve./minLuminanceBanksCurve;
    banksSquareRootLaw2 = maxLuminanceBanksCurve./midLuminanceBanksCurve;
    
    for lumIndex = 1:nLuminances
        cpd = figData.banksEtAlReplicate.cyclesPerDegree(lumIndex,:);
        thresholdContrasts = [figData.banksEtAlReplicate.mlptThresholds(lumIndex,:).thresholdContrasts];
        contrastSensitivities(lumIndex,:) = 1./(thresholdContrasts*referenceContrast);
        plot(cpd,contrastSensitivities(lumIndex,:), 'ko', ...
            'MarkerEdgeColor', squeeze(colors(lumIndex,:)), ...
            'MarkerFaceColor', min(1, squeeze(colors(lumIndex,:)) + [0.6 0.6 0.6]), ...
            'MarkerSize',11,'LineWidth',2);  
    end
    
    dataSquareRootLaw1 = squeeze(contrastSensitivities(2,:)) ./ squeeze(contrastSensitivities(1,:))
    dataSquareRootLaw2 = squeeze(contrastSensitivities(3,:)) ./ squeeze(contrastSensitivities(2,:))
    
    hL = legend({sprintf('  %2.1f cd/m2', luminances(1)), sprintf(' %2.1f cd/m2', luminances(2)), sprintf('%2.1f cd/m2', luminances(3))});
    set(hL, 'FontName', 'Monaco');
    set(gca, 'YScale', 'log', 'XScale', 'log', 'YTick', csTicks, 'XTick', sfTicks, ...
        'XLim', [sfLims(1)-dx1 sfLims(2)+dx2], 'YLim', [csLims(1)-dy1 csLims(2)+dy2]);
    box on; grid on;
    set(gca, 'FontSize', 14, 'TickLength',[0.02, 0.02]);
    xlabel('spatial frequency (c/deg)', 'FontWeight', 'bold');
    ylabel('contrast sensitivity', 'FontWeight', 'bold');
    set(gca, 'Color', 'none');
    set(hFig, 'Color', 'none');
    
    NicePlot.exportFigToPDF('GeislerOptics2mm.pdf', hFig, 300);
    
    yLim =   [2.0 5];
    yTicks = 1:0.5:5;
    
    hFig = figure(2);
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 400 400]);
    clf;
    hold on;
    plot(cpd, dataSquareRootLaw1, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0.75 0.75 0.0], 'LineWidth', 1.5);
    plot(cpd, dataSquareRootLaw2, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0.0 0.75 0.75], 'LineWidth', 1.5);
    plot([1 200], sqrt(10)+[0 0], 'k-', 'LineWidth', 1.5);
    plot(hiResSF,banksSquareRootLaw1, '-', 'Color', [0.75 0.75 0.0], 'LineWidth', 2);
    plot(hiResSF,banksSquareRootLaw2, '-' ,'Color', [0.0 0.75 0.75], 'LineWidth', 2);
    
    hL = legend({sprintf(' %2.1f :  %-2.1f cd/m2', luminances(2),luminances(1)), sprintf('%-2.1f : %-2.1f cd/m2', luminances(3), luminances(2)), 'sq.root(10)'})
    set(hL, 'FontName', 'Monaco');
    set(gca, 'YScale', 'log', 'XScale', 'log', 'YTick', yTicks', 'YTickLabels', sprintf('%2.2f\n',yTicks), 'XTick', sfTicks, ...
        'XLim', [sfLims(1)-dx1 sfLims(2)+dx2], 'YLim', yLim);
    box on; grid on;
    set(gca, 'FontSize', 14, 'TickLength',[0.02, 0.02]);
    xlabel('spatial frequency (c/deg)', 'FontWeight', 'bold');
    ylabel('cs ratio', 'FontWeight', 'bold');
    NicePlot.exportFigToPDF('GeislerOptics2mmRatios.pdf', hFig, 300);
    
end