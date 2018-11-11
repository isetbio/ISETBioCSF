function visualizeResponses(coneMosaic, coneExcitations, figNo)
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1420 740], 'Color', [1 1 1], 'Name', 'Mosaic and 3 response instances');
    axHandle = subplot(2,2,1);
    coneMosaic.visualizeGrid('axesHandle', axHandle, 'displayVisualDegs', true);
    xlabel('space (degs', 'FontWeight', 'bold');
    ylabel('space (degs', 'FontWeight', 'bold');
    title('cone mosaic');
    
    for k = 1:3
        axHandle = subplot(2,2,1+k);
        coneMosaic.renderActivationMap(axHandle, squeeze(coneExcitations(k,:,:)), ...
                'mapType', 'modulated disks', ...
                'signalRange', [0 max(coneExcitations(:))], ...
                'showColorBar', ~true, ...
                'labelColorBarTicks', ~true);
        set(gca, 'XTickLabel', {}, 'YTickLabel', {});
        xlabel(''); ylabel('');
        title(sprintf('response instance %d', k));
    end
end