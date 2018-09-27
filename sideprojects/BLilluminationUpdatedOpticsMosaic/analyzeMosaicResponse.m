function analyzeMosaicResponse
    load('data.mat', 'theHexMosaic', 'theScene', 'theOI', 'absorptions');
    
    theHexMosaic
    pause
    theHexMosaic.pigment
    pause
    size(theHexMosaic.pattern)
    size(absorptions)
    idx = find(theHexMosaic.pattern(:) == 2);
    absorptions = absorptions(:);
    LconeMeanActivation = absorptions(idx);
    
    idx = find(theHexMosaic.pattern(:) == 3);
    absorptions = absorptions(:);
    MconeMeanActivation = absorptions(idx);
    
    idx = find(theHexMosaic.pattern(:) == 4);
    absorptions = absorptions(:);
    SconeMeanActivation = absorptions(idx);
    
    figure(1); clf;
    subplot(1,3,1);
    histogram(LconeMeanActivation, 0:2:100, 'FaceColor', [0.9 0.9 0.9]);
    hold on;
    plot(mean(LconeMeanActivation)*[1 1], [0 1000], 'r-', 'LineWidth', 1.5);
    set(gca, 'XLim', [0 60], 'YLim', [0 600], 'FontSize', 12);
    title(sprintf('mean L cone excitations: %2.2f R*/%2.0f ms bin', mean(LconeMeanActivation), theHexMosaic.integrationTime*1000));
    xlabel(sprintf('excitations (R*/%2.0f ms bin)', theHexMosaic.integrationTime*1000), 'FontWeight', 'Bold');
    
    subplot(1,3,2);
    histogram(MconeMeanActivation, 0:2:100, 'FaceColor', [0.9 0.9 0.9]);
    hold on;
    plot(mean(MconeMeanActivation)*[1 1], [0 1000], 'r-', 'LineWidth', 1.5);
    set(gca, 'XLim', [0 60], 'YLim', [0 600], 'FontSize', 12);
    title(sprintf('mean M cone excitations: %2.2f R*/%2.0f ms bin', mean(MconeMeanActivation), theHexMosaic.integrationTime*1000));
    xlabel(sprintf('excitations (R*/%2.0f ms bin)', theHexMosaic.integrationTime*1000), 'FontWeight', 'Bold');
    
    subplot(1,3,3);
    histogram(SconeMeanActivation, 0:2:60, 'FaceColor', [0.9 0.9 0.9]);
    hold on;
    plot(mean(SconeMeanActivation)*[1 1], [0 1000], 'r-', 'LineWidth', 1.5);
    set(gca, 'XLim', [0 60], 'YLim', [0 600], 'FontSize', 12);
    title(sprintf('mean S cone excitations: %2.2f R*/%2.0f ms bin', mean(SconeMeanActivation), theHexMosaic.integrationTime*1000));
    xlabel(sprintf('excitations (R*/%2.0f ms bin)', theHexMosaic.integrationTime*1000), 'FontWeight', 'Bold');
    
end

