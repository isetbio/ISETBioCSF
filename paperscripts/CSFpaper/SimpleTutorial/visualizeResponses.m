function visualizeResponses(coneMosaic, coneExcitations, figNo)
    figure(figNo); clf;
    axHandle = subplot(2,2,1);
    coneMosaic.visualizeGrid('axesHandle', axHandle, 'displayVisualDegs', true);
    
    for k = 1:3
        axHandle = subplot(2,2,1+k);
        coneMosaic.renderActivationMap(axHandle, squeeze(coneExcitations(k,:,:)), 'mapType', 'modulated disks');
    end
end