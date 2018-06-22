function [noStimData, stimData] = transformDataWithSeparableSpaceTimeComponents(noStimData, stimData, thresholdParams, paramsList, visualizeSignals)
% [noStimData, stimData] = transformDataWithSeparableSpaceTimeComponents(noStimData, stimData, thresholdParams)
%

fprintf('Transforming data via projection to %d spatial components (PCA)\n', thresholdParams.PCAComponents);

repsDimension = 1;
spatialDimension = 2;
temporalDimension = 3;

% Subtract the noStimData so that zero modulation gives zero response for both isomerizations and photocurrrents
[noStimData, stimData] = subtractMeanOfNoStimData(noStimData, stimData, thresholdParams.signalSource, repsDimension, temporalDimension);

if (strcmp(thresholdParams.signalSource,'photocurrents'))
    if (ndims(noStimData.responseInstanceArray.theMosaicPhotocurrents) ~= 3)
        error('transformDataWithSeparableSpaceTimeComponents not yet implemented for other than 3D response arrays, ndims = %d\n', ndims(noStimData.responseInstanceArray.theMosaicPhotocurrents));
    end
    repsNum = size(noStimData.responseInstanceArray.theMosaicPhotocurrents,repsDimension);
	d = noStimData.responseInstanceArray.theMosaicPhotocurrents;
    d = cat(repsDimension, d, stimData.responseInstanceArray.theMosaicPhotocurrents);
    noStimData.responseInstanceArray.theMosaicPhotocurrents = [];
    stimData.responseInstanceArray.theMosaicPhotocurrents = [];
else
    if (ndims(noStimData.responseInstanceArray.theMosaicIsomerizations) ~= 3)
        error('transformDataWithSeparableSpaceTimeComponents not yet implemented for other than 3D response arrays, ndims = %d\n', ndims(noStimData.responseInstanceArray.theMosaicIsomerizations));
    end
    repsNum = size(noStimData.responseInstanceArray.theMosaicIsomerizations,repsDimension);
	d = noStimData.responseInstanceArray.theMosaicIsomerizations;
	d = cat(repsDimension, d, stimData.responseInstanceArray.theMosaicIsomerizations);
    noStimData.responseInstanceArray.theMosaicIsomerizations = [];
    stimData.responseInstanceArray.theMosaicIsomerizations = [];
end

spatialBinsNum = size(d, spatialDimension);
timeBinsNum = size(d, temporalDimension);

% Find PCAs for the ensemble of STIM and NOSTIM data
d = reshape(permute(d, [1 3 2]), [repsNum*2*timeBinsNum spatialBinsNum]);
[d,~,varianceExplained] = transformDataWithPCA(d, thresholdParams.PCAComponents, thresholdParams.STANDARDIZE);
for compIndex = 1:thresholdParams.PCAComponents
   fprintf('Variance explained by component #%d: %2.2f (accum: %2.2f)\n', compIndex, varianceExplained(compIndex), sum(varianceExplained(1:compIndex)));
end
        
d = permute(reshape(d, [repsNum*2 timeBinsNum thresholdParams.PCAComponents]), [1 3 2]);

if (strcmp(thresholdParams.signalSource,'photocurrents'))
	noStimData.responseInstanceArray.theMosaicPhotocurrents = d(1:repsNum,:,:);
	stimData.responseInstanceArray.theMosaicPhotocurrents = d(repsNum+1:end,:,:);

    photocurrentsBasedPCAresponseRange = [...
        min([min(stimData.responseInstanceArray.theMosaicPhotocurrents(:)) min(noStimData.responseInstanceArray.theMosaicPhotocurrents(:))]) ...
        max([max(stimData.responseInstanceArray.theMosaicPhotocurrents(:)) max(noStimData.responseInstanceArray.theMosaicPhotocurrents(:))])];
else
	noStimData.responseInstanceArray.theMosaicIsomerizations = d(1:repsNum,:,:);
	stimData.responseInstanceArray.theMosaicIsomerizations = d(repsNum+1:end,:,:);
    isomerizationsBasedPCAresponseRange = [ ...
        min([min(stimData.responseInstanceArray.theMosaicIsomerizations(:)) min(noStimData.responseInstanceArray.theMosaicIsomerizations(:))]) 
        max([max(stimData.responseInstanceArray.theMosaicIsomerizations(:)) max(noStimData.responseInstanceArray.theMosaicIsomerizations(:))])];
end


if (visualizeSignals) 
    visualizedPCAcomponentsNum = min([thresholdParams.PCAComponents 3]);
    hFig = figure(1234); clf;
    set(hFig, 'Position', [10 10 400*visualizedPCAcomponentsNum 800]);
    
    for pcaComponentIndex = 1:visualizedPCAcomponentsNum
        
        subplot(2, visualizedPCAcomponentsNum, pcaComponentIndex);
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            plot(noStimData.responseInstanceArray.timeAxis, squeeze(noStimData.responseInstanceArray.theMosaicPhotocurrents(:,pcaComponentIndex,:)), 'k-');
            set(gca, 'YLim', photocurrentsBasedPCAresponseRange,  'XLim', [noStimData.responseInstanceArray.timeAxis(1) noStimData.responseInstanceArray.timeAxis(timeBinsNum)]);
            title(sprintf('NO-STIM\nphotocurrents-based component-%d response', pcaComponentIndex));
        else
            plot(noStimData.responseInstanceArray.timeAxis, squeeze(noStimData.responseInstanceArray.theMosaicIsomerizations(:,pcaComponentIndex,:)), 'k-');
            set(gca, 'YLim', isomerizationsBasedPCAresponseRange,  'XLim', [noStimData.responseInstanceArray.timeAxis(1) noStimData.responseInstanceArray.timeAxis(timeBinsNum)]);
            title(sprintf('NO-STIM\nisomerizations-based component-%d response', pcaComponentIndex));
        end
        if (pcaComponentIndex == 1)
            ylabel('V1 filter bank energy');
        end
        set(gca, 'FontSize', 14);
        
        subplot(2, visualizedPCAcomponentsNum, pcaComponentIndex+visualizedPCAcomponentsNum);
        if (strcmp(thresholdParams.signalSource,'photocurrents'))
            plot(stimData.responseInstanceArray.timeAxis, squeeze(stimData.responseInstanceArray.theMosaicPhotocurrents(:,pcaComponentIndex,:)), 'k-');
            set(gca, 'YLim', photocurrentsBasedPCAresponseRange,  'XLim', [noStimData.responseInstanceArray.timeAxis(1) noStimData.responseInstanceArray.timeAxis(timeBinsNum)]);
            title(sprintf('C = %2.5f%%\nphotocurrents-based component-%d response', stimData.testContrast*100, pcaComponentIndex));
        else
            plot(stimData.responseInstanceArray.timeAxis, squeeze(stimData.responseInstanceArray.theMosaicIsomerizations(:,pcaComponentIndex,:)), 'k-');
            set(gca, 'YLim', isomerizationsBasedPCAresponseRange,  'XLim', [noStimData.responseInstanceArray.timeAxis(1) noStimData.responseInstanceArray.timeAxis(timeBinsNum)]);
            title(sprintf('C = %2.5%%\nisomerizations-based component-%d response', stimData.testContrast*100, pcaComponentIndex));
        end
        if (pcaComponentIndex == 1)
            ylabel('V1 filter bank energy');
        end
        xlabel('time (ms)'); 
        set(gca, 'FontSize', 14);
        
        drawnow;  
    end % for pcaComponentIndex 
    
    % Save figure
    theProgram = mfilename;
    rwObject = IBIOColorDetectReadWriteBasic;
    data = 0;
    fileName = sprintf('%s-based_%s_outputs', thresholdParams.signalSource, thresholdParams.method);
    rwObject.write(fileName, data, paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');

end % visualize transformed signals

end