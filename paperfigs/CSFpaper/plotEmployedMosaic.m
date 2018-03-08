function plotEmployedMosaic()
    %load('/Users/nicolas/Desktop/coneMosaic50cpd.mat');
    %theConeMosaic = theData;
    
    load('/Users/nicolas/Desktop/coneMosaicEccLMS3.5degs.mat');
    theConeMosaic = theData;
    xLimsDegs = 1.75*[-1 1];
    yLimsDegs = 1.75*[-1 1];
    
%     load('/Users/nicolas/Desktop/coneMosaicEccLMS0.98degs.mat');
%     theConeMosaic = theData;
%     xLimsDegs = 0.5*[-1 1];
%     yLimsDegs = 0.5*[-1 1];
    
    %load('theHexMosaic0.75degs.mat');
    %theConeMosaic = theHexMosaic;
    %xLimsDegs = 0.35*[-1 1];
    %yLimsDegs = 0.18*[-1 1];
    
    hFig = figure(1);
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 940 1025]);
    ax = subplot('Position', [0.03 0.05 0.96 0.93]);
    theConeMosaic.visualizeGrid('axesHandle', ax, 'apertureShape', 'disks', ...
         'ForegroundColor', 'k',...
         'overlayHexMesh', false, ...
         'overlayConeDensityContour', 'theoretical_and_measured', ...
         'coneDensityContourLevels', 1e3*[35 42 50 60 80 110 150 220], ...
         'labelConeTypes', false);
    
%     theConeMosaic.visualizeGrid('axesHandle', gca, 'apertureShape', 'disks', ...
%          'ForegroundColor', 'k',...
%          'overlayHexMesh', true, ...
%          'labelConeTypes', true);
     
    
    xLimsMeters = xLimsDegs*300 *1e-6;
    yLimsMeters = yLimsDegs*300 *1e-6;
    
    xTicksDegs = 0.5*[-10:1:10];
    xTicks = xTicksDegs*300 *1e-6;
    xTickLabels = sprintf('%-2.2f\n', xTicksDegs);
    box on;
    set(gca, 'FontSize', 14, 'XLim', xLimsMeters, ...
        'YLim', yLimsMeters, ...
        'YTick', xTicks, 'YTickLabel', {},...
        'XTick', xTicks, 'XTickLabel', xTickLabels);
    NicePlot.exportFigToPDF('EccLMSMosaic.pdf', hFig, 300);
    xlabel('space (degs)');
end