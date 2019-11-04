function generateFigureDorDavidsCappadociaTalk
    generateFig1()
    generateFig2()
end


function generateFig2()
    load('Fig2')
    expData.x = Fig2(:,1)
    expData.y = Fig2(:,2)
    expData.errorBar = (Fig2(:,4)-expData.y)*0.7;
    poolingData.x = Fig2(:,5);
    poolingData.y = Fig2(:,6);
    tmp= [...
-8.280802292264E-1	6.724137931034E-1;
-7.550143266476E-1	6.853448275862E-1;
-6.819484240688E-1	6.939655172414E-1;
-6.002865329513E-1	7.068965517241E-1;
-5.100286532951E-1	7.241379310345E-1;
-4.197707736390E-1	7.413793103448E-1;
-3.295128939828E-1	7.586206896552E-1;
-2.177650429799E-1	7.758620689655E-1;
-1.404011461318E-1	7.974137931034E-1;
-4.154727793696E-2	8.275862068966E-1;
3.151862464183E-2	8.448275862069E-1;
1.174785100287E-1	8.663793103448E-1;
1.905444126074E-1	8.836206896552E-1;
2.979942693410E-1	9.137931034483E-1;
3.968481375358E-1	9.568965517241E-1;
5.000000000000E-1	1.004310344828E0;
5.730659025788E-1	1.025862068966E0;
6.418338108883E-1	1.060344827586E0;
7.277936962751E-1	1.073275862069E0;
7.879656160458E-1	1.099137931034E0;
8.481375358166E-1	1.125000000000E0;
9.212034383954E-1	1.159482758621E0;
9.985673352436E-1	1.185344827586E0;
1.088825214900E0	1.219827586207E0;
1.183381088825E0	1.262931034483E0;
1.269340974212E0	1.297413793103E0;
1.342406876791E0	1.336206896552E0;
1.419770773639E0	1.375000000000E0;
1.492836676218E0	1.413793103448E0;
1.570200573066E0	1.452586206897E0;
1.660458452722E0	1.495689655172E0;
1.737822349570E0	1.521551724138E0;
1.819484240688E0	1.560344827586E0;];
    noSummationData.x = tmp(:,1);
    noSummationData.y = tmp(:,2);

    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 845 425], 'Color', [1 1 1]);
    errorbar(expData.x, expData.y, expData.errorBar, 'k-', 'MarkerSize', 1, 'LineWidth', 3.0);
    hold on;
    p1 = plot(expData.x, expData.y, 'ko-', 'MarkerSize', 16, 'MarkerFaceColor', [0.3 0.3 0.3], 'LineWidth', 3.0); 
    p2 = plot(poolingData.x, poolingData.y, 'r:', 'LineWidth', 4.0);
    p3 = plot(noSummationData.x, noSummationData.y, 'b-', 'LineWidth', 4.0, 'Color', [0.5 0.5 0.5]);
    
    legend([p1 p2 p3], {'mean data \pm 2 SD (condition 3)', 'Gaussian \sigma: 1.7 arcmin', 'no summation'}, 'Location', 'NorthWest')
    set(gca, 'XLim', [-1 2], 'YLim', [0.5 2], 'LineWidth', 1.0);
    set(gca, 'XTick', -1:1:2, 'YTick', 0:1:2);
    set(gca, 'FontSize', 30);
    grid on; box on
    xlabel('log stimulus area (arcmin^2)');
    ylabel('threshold energy');
    
    drawnow
    NicePlot.exportFigToPDF('fig2.pdf', hFig, 300);
end


function generateFig1()
    load('Fig1')
    expData.x = Fig1(:,1);
    expData.y = Fig1(:,2);
    expData.errorBar = (Fig1(:,8)-expData.y)*0.7;
    poolingData.x = Fig1(:,3);
    poolingData.y = Fig1(:,4);
    noSummationData.x = Fig1(:,5);
    noSummationData.y = Fig1(:,6);
    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 845 425], 'Color', [1 1 1]);
    errorbar(expData.x, expData.y, expData.errorBar, 'k-', 'MarkerSize', 1,  'LineWidth', 3.0); 
    hold on;
    p1 = plot(expData.x, expData.y, 'ko-', 'MarkerSize', 16, 'MarkerFaceColor', [0.3 0.3 0.3], 'LineWidth', 3.0); 
    p2 = plot(poolingData.x, poolingData.y, 'r:', 'LineWidth', 4.0);
    p3 = plot(noSummationData.x, noSummationData.y, 'b-', 'LineWidth', 4.0, 'Color', [0.5 0.5 0.5]);
    legend([p1 p2 p3], {'mean data \pm 2 SD (condition 1)', 'Gaussian \sigma: 1.7 arcmin', 'no summation'}, 'Location', 'NorthWest')
    set(gca, 'XLim', [-1 2], 'YLim', [0.5 2], 'LineWidth', 1.0);
    set(gca, 'XTick', -1:1:2, 'YTick', 0:1:2);
    set(gca, 'FontSize', 30);
    grid on; box on
    xlabel('log stimulus area (arcmin^2)');
    ylabel('threshold energy');
    
    drawnow
    NicePlot.exportFigToPDF('fig1.pdf', hFig, 300);
end

