function visualizeAllResponses(stimDescriptor, theMosaic, contrastLevel, theInstance, nTrials, trialNo, resourcesDir, figNo)
    
    if (strcmp(stimDescriptor, 'zeroContrast'))
        fname = fullfile(resourcesDir, sprintf('zeroContrast_nTrials_%d.mat',  nTrials));
    else
        fname = fullfile(resourcesDir, sprintf('%s_contrast_%2.4f_instance_%1.0f_nTrials_%d.mat', stimDescriptor, contrastLevel, theInstance, nTrials));
    end
    load(fname, 'coneExcitations', 'photoCurrents', 'eyeMovementPaths',  'emPathsDegs', 'timeAxis')

    visualizeDynamicResponse(theMosaic, coneExcitations, photoCurrents, eyeMovementPaths, emPathsDegs, timeAxis, stimDescriptor, trialNo, figNo);    
end

function visualizeDynamicResponse(theMosaic, coneExcitations, photoCurrents, eyeMovementPaths, emPathsDegs, timeAxis, stimDescriptor, trialNo, figNo)
    
    [idxOfConesAlongHorizMeridian, idxOfConesAlongVertMeridian, ...
            eccDegsOfConesAlongHorizMeridian, ...
            eccDegsOfConesAlongVertMeridian, ...
            idxOfConesAlongHorizMeridianInSerializedList, idxOfConesAlongVertMeridianInSerializedList, ] = theMosaic.indicesForConesAlongMeridians();
        

    singleTrialConeExcitationResponse = squeeze(coneExcitations(trialNo,:,:));
    coneExcitationResponseRange = prctile(singleTrialConeExcitationResponse(:), [5 95]);
    
    singleTrialPhotocurrentResponse = squeeze(photoCurrents(trialNo,:,:));
    photocurrentResponseRange = prctile(singleTrialPhotocurrentResponse(:), [5 95]);
    
    singleTrialEyeMovement = squeeze(eyeMovementPaths(trialNo,:,:));
    singleTrialEyeMovementDegs = squeeze(emPathsDegs(trialNo,:,:));
    
    hFig = figure(figNo);
    
    for timeBin = 1:size(singleTrialConeExcitationResponse,2)
        axHandle = subplot(2,5,6);
        theMosaic.renderActivationMap(axHandle, squeeze(singleTrialConeExcitationResponse(:,timeBin)), ...
                'mapType', 'modulated disks', ...
                'signalRange', coneExcitationResponseRange, ...
                'showColorBar', true, ...
                'labelColorBarTicks', true, ...
                'titleForColorBar', 'R*/cone/tau');

        % The eye movement trajectory
        subplot(2,5,1);
        plot(singleTrialEyeMovementDegs(1:timeBin,1),  singleTrialEyeMovementDegs(1:timeBin,2), 'k-'); hold on;
        plot(singleTrialEyeMovementDegs(timeBin,1)+ [-0.1 0.1],  singleTrialEyeMovementDegs(timeBin,2)*[1 1], 'b-', 'LineWidth', 1.5);
        plot(singleTrialEyeMovementDegs(timeBin,1)*[1 1],  singleTrialEyeMovementDegs(timeBin,2)+ [-0.1 0.1], 'b-', 'LineWidth', 1.5); 
        hold off;
        axis 'square';
        set(gca, 'YLim', [-0.2 0.2], 'XLim', [-0.2 0.2]);
        xlabel('horizontal position (degs)');
        ylabel('vertical position (degs)');
        
        % The cone excitations for cones along the vertical meridian
        subplot(2,5,2);
        responseVector = squeeze(singleTrialConeExcitationResponse(idxOfConesAlongVertMeridianInSerializedList,timeBin));
        plot(eccDegsOfConesAlongVertMeridian, responseVector, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8]);
        axis 'square';
        set(gca, 'YLim', coneExcitationResponseRange);
        xlabel('vertical position (degs)');
        ylabel('cone excitations (R*/cone/tau)');
        title('cone excitation responses (vert meridian)');
        
         % The cone excitations for cones along the horizontal meridian
        subplot(2,5,3);
        responseVector = squeeze(singleTrialConeExcitationResponse(idxOfConesAlongHorizMeridianInSerializedList,timeBin));
        plot(eccDegsOfConesAlongHorizMeridian, responseVector, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8]);
        axis 'square';
        set(gca, 'YLim', coneExcitationResponseRange);
        xlabel('horizontal position (degs)');
        ylabel('cone excitations (R*/cone/tau)');
        title('cone excitation responses (horiz meridian)');
        
        % XT plot of cone excitations along the horizontal meridian
        subplot(2,5,4);
        imagesc(eccDegsOfConesAlongHorizMeridian, timeAxis*1000, (singleTrialConeExcitationResponse(idxOfConesAlongHorizMeridianInSerializedList,:))');
        hold on;
         % superimpose x-eye movement trajectory
        plot(singleTrialEyeMovementDegs(:,1), timeAxis*1000, 'r-', 'LineWidth', 1.5);
        set(gca, 'YTick', 0:100:1000);
        hold off
        axis 'square';
        axis 'xy'
        xlabel('horizontal position (degs)');
        ylabel('time (msec)');
        colormap(gray);
        title('cone excitations (horiz meridian)');
        
        % XT plot of cone excitations along the vertical meridian
        subplot(2,5,5);
        imagesc(eccDegsOfConesAlongVertMeridian, timeAxis*1000, (singleTrialConeExcitationResponse(idxOfConesAlongVertMeridianInSerializedList,:))');
        hold on;
        % superimpose y-eye movement trajectory
        plot(-singleTrialEyeMovementDegs(:,2), timeAxis*1000, 'g-', 'LineWidth', 1.5);
        set(gca, 'YTick', 0:100:1000);
        hold off
        axis 'square';
        axis 'xy'
        xlabel('vertical position (degs)');
        ylabel('time (msec)');
        colormap(gray);
        title('cone excitations (vert meridian)');
        
        
        % The photocurrents for cones along the vertical meridian
        subplot(2,5,7);
        responseVector = squeeze(singleTrialPhotocurrentResponse(idxOfConesAlongVertMeridianInSerializedList,timeBin));
        plot(eccDegsOfConesAlongVertMeridian, responseVector, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8]);
        axis 'square';
        set(gca, 'YLim', photocurrentResponseRange);
        xlabel('vertical position (degs)');
        ylabel('photocurrent (pAmps)');
        title('photocurrent responses (vert meridian)');
        
        % The photocurrents for cones along the horizontal meridian
        subplot(2,5,8);
        responseVector = squeeze(singleTrialPhotocurrentResponse(idxOfConesAlongHorizMeridianInSerializedList,timeBin));
        plot(eccDegsOfConesAlongHorizMeridian, responseVector, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8]);
        axis 'square';
        set(gca, 'YLim', photocurrentResponseRange);
        xlabel('horizontal position (degs)');
        ylabel('photocurrent (pAmps)');
        title('photocurrent responses (horiz meridian)');
        
        % XT plot of photocurrents along the horizontal meridian
        subplot(2,5,9);
        imagesc(eccDegsOfConesAlongHorizMeridian, timeAxis*1000, (singleTrialPhotocurrentResponse(idxOfConesAlongHorizMeridianInSerializedList,:))');
        hold on;
        % superimpose x-eye movement trajectory
        plot(singleTrialEyeMovementDegs(:,1), timeAxis*1000, 'r-', 'LineWidth', 1.5);
        set(gca, 'YTick', 0:100:1000);
        hold off
        axis 'square';
        axis 'xy'
        xlabel('horizontal position (degs)');
        ylabel('time (msec)');
        title('photocurrent (horiz meridian)');
        
        % XT plot of photocurrents along the vertical meridian
        subplot(2,5,10);
        imagesc(eccDegsOfConesAlongVertMeridian, timeAxis*1000, (singleTrialPhotocurrentResponse(idxOfConesAlongVertMeridianInSerializedList,:))');
        hold on;
        % superimpose y-eye movement trajectory
        plot(-singleTrialEyeMovementDegs(:,2), timeAxis*1000, 'g-', 'LineWidth', 1.5);
        set(gca, 'YTick', 0:100:1000);
        hold off
        axis 'square';
        axis 'xy'
        xlabel('vertical position (degs)');
        ylabel('time (msec)');
        title('photocurrent (vert meridian)');
        
        colormap(gray);
        
        drawnow;
    end
end

