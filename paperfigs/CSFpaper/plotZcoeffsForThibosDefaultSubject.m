function plotZcoeffsForThibosDefaultSubject
Zcoeffs = [
   0.194587472215542
   -0.013733074853016
   0.144577423523313
  -0.038004150615687
  -0.004193502228939
  -0.047213114954317
  -0.018800220466921
  -0.005811359533760
   0.001381338758501
   0.000613179115350
  -0.003901639791644
   0.003668859399324
   0.010573866630993
  -0.003212468492853
  -0.003862587927572
  ];
    
    
    ZcoeffNames = {...
'piston'
'vertical_tilt'
'horizontal_tilt'
'oblique_astigmatism'
'defocus'
'vertical_astigmatism'
'vertical_trefoil'
'vertical_coma'
'horizontal_coma'
'oblique_trefoil'
'oblique_quadrafoil'
'oblique_secondary_astigmatism'
	'spherical'
	'vertical_secondary_astigmatism'
	'vertical_quadrafoil'
};

for k = 1:numel(ZcoeffNames)
    ZcoeffNames2{k} = sprintf('%s (%d)', strrep(ZcoeffNames{k}, '_', ' '), k);
end
[radialOrder, meridionalFrequency] = wvfOSAIndexToZernikeNM(1:15)
pause

    hFig = figure(1);clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 620 620]);
    %subplot(1,2,1);
    stem(1:numel(Zcoeffs), Zcoeffs, 'r', 'filled', 'MarkerSize', 14, 'LineWidth', 1.5);
    set(gca, 'XLim', [0 16], 'YLim', 0.1*[-2.5 2.5], 'XTick', 1:1:16, 'YTick', -3:0.1:3);
    grid on; box on
    set(gca, 'XTick', 1:15, 'XTickLabel', ZcoeffNames2, 'FontSize', 14);
    xlabel('coefficient', 'FontWeight', 'bold', 'FontSize', 16);
    ylabel('c_k', 'FontWeight', 'bold', 'FontSize', 16);
    xtickangle(90);
    NicePlot.exportFigToPDF('ZernikeCoeffs.pdf', hFig, 300);
    
%     subplot(1,2,2)
%     map2D = zeros(6,11);
%     for k = 1:numel(radialOrder)
%         row = 1+radialOrder(k);
%         col = meridionalFrequency(k) + 6;
%         map2D(row,col) = Zcoeffs(k);
%     end
%     radialOrderAxis = 1:6;
%     meridialFrequencyAxis = -5:5;
%     imagesc(meridialFrequencyAxis, radialOrderAxis, map2D);
%     colormap(bone(256));
%     set(gca, 'CLim', 0.1*[-2.0 2.0]);
%     axis 'image'
 
    
end

