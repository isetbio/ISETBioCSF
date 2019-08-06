function visualizeResponses(coneMosaic, coneExcitations, figNo)
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1500 950], 'Color', [1 1 1], 'Name', 'Mosaic and 3 response instances');
    axHandle = subplot(2,3,2);
    coneMosaic.visualizeGrid('axesHandle', axHandle, 'ticksInVisualDegs', true);
    xlabel('space (degs)', 'FontWeight', 'bold');
    ylabel('space (degs)', 'FontWeight', 'bold');
    title('cone mosaic');
    
    for k = 1:3
        axHandle = subplot(2,3,3+k);
        if (k == 3)
            coneMosaic.renderActivationMap(axHandle, squeeze(coneExcitations(k,:,:)), ...
                'mapType', 'modulated disks', ...
                'signalRange', [0 max(coneExcitations(:))], ...
                'showColorBar', true, ...
                'labelColorBarTicks', true, ...
                'titleForColorBar', sprintf('R*/cone/%2.0fms', coneMosaic.integrationTime*1000));
        else 
            coneMosaic.renderActivationMap(axHandle, squeeze(coneExcitations(k,:,:)), ...
                'mapType', 'modulated disks', ...
                'signalRange', [0 max(coneExcitations(:))], ...
                'showColorBar', ~true, ...
                'labelColorBarTicks', ~true);
        end
        if (k > 1)
        	set(gca, 'YTickLabel', {});
            ylabel('');
        end
        xlabel('space (degs)', 'FontWeight', 'bold'); 
        title(sprintf('response instance #%d', k));
    end
end