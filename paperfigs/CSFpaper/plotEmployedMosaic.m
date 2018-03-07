function plotEmployedMosaic()
    %load('/Users/nicolas/Desktop/coneMosaic50cpd.mat');
    %theConeMosaic = theData;
    
    load('theHexMosaic0.75degs.mat');
    theConeMosaic = theHexMosaic;
    
    hFig = figure();
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1200 600]);
    theConeMosaic.visualizeGrid('axesHandle', gca, 'apertureShape', 'disks', ...
         'ForegroundColor', 'k',...
        'overlayHexMesh', false, ...
        'labelConeTypes', true);
    
    xLimsDegs = 0.35*[-1 1];
    xLimsMeters = xLimsDegs*300 *1e-6;
    yLimsDegs = 0.18*[-1 1];
    yLimsMeters = yLimsDegs*300 *1e-6;
    
    xTicksDegs = 0.1*[-10:1:10];
    xTicks = xTicksDegs*300 *1e-6;
    xTickLabels = sprintf('%-2.2f\n', xTicksDegs);
    box on;
    set(gca, 'FontSize', 14, 'XLim', xLimsMeters, ...
        'YLim', yLimsMeters, ...
        'YTick', xTicks, 'YTickLabel', {},...
        'XTick', xTicks, 'XTickLabel', xTickLabels);
    NicePlot.exportFigToPDF('HexMosaicNoMesh.pdf', hFig, 300);
    xlabel('space (degs)');
end