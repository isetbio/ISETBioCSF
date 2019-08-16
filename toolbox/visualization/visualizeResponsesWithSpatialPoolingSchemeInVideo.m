function visualizeResponsesWithSpatialPoolingSchemeInVideo(filterBankName, filterBank, colorModulationParams, spatialParams, mosaicParams, topLevelDirParams, thresholdParams, constantParamsList)
    
    % Make sure we have a valid filterBankName
    if (~ismember(filterBankName, {'svmV1FilterBank', 'svmV1FilterEnsemble'}))
        error('filterBankName must be either ''svmV1FilterBank'' *OR* ''svmV1FilterEnsemble''. It is: ''%s''.', filterBankName);
    end

    % Load the employed mosaic
    theMosaic = loadTheMosaic(topLevelDirParams, mosaicParams);
    
    % Load response data
    [stimData, noStimData] = loadResponseData(colorModulationParams, constantParamsList);
    
    % Extract isomerizations data for video
    signalSource = 'isomerizations'; 
    videoData = extractVideoData(stimData, noStimData, filterBank, theMosaic, spatialParams, thresholdParams, signalSource);
    
    save('tmp.mat', 'theMosaic', 'videoData');
    fprintf('saved data');
    pause
    % Render video for isomerizations
    videoFileName = sprintf('%sVideo', signalSource);
    generateVideo(videoData, theMosaic, videoFileName);
    
    
    % Extract photocurrent data for video
    signalSource = 'photocurrents';
    videoData = extractVideoData(stimData, noStimData, filterBank, theMosaic, spatialParams, thresholdParams, signalSource);
    
    % Render video for photocurrents
    videoFileName = sprintf('%sVideo', signalSource);
    generateVideo(videoData, theMosaic, videoFileName);

end

function generateVideo(videoData, theMosaic, videoFileName)
    
    % only visualize from 0 to 100 msec
    visualizedTimeRangeMilliseconds = [0 125];
    tBinIndices = find(videoData.nullResponseInstances.timeAxis>=visualizedTimeRangeMilliseconds(1)/1000 & videoData.nullResponseInstances.timeAxis <= visualizedTimeRangeMilliseconds(2)/1000);
    
    % only visualize portion of the mosaic
    visualizedMosaicRangeDegs = 0.5*[-1 1];
    visualizedMosaicRange = 0.5*visualizedMosaicRangeDegs * theMosaic.micronsPerDegree * 1e-6;
    
    % visualized response range
    responseRange = prctile(videoData.differentialResponses(:), 50+45*[-1.0 1.0]);
    
    % trials number
    nTrials = size(videoData.nullResponseInstances.theMosaicIsomerizations,1);
    
    % Gamma corrected gray colormap
    responseColorMap = brewermap(1024, '*Greys');
    responseColorMap = responseColorMap .^ (1/2.2);
    
    % Open video stream
    videoOBJ = VideoWriter(videoFileName, 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 20;
    videoOBJ.Quality = 100;
    videoOBJ.open();

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1230 410], 'Color', [1 1 1]);
    
    % Subfig layout
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', 3, ...
       'heightMargin', 0.08, ...
       'widthMargin', 0.08, ...
       'leftMargin', 0.06, ...
       'rightMargin', 0.02, ...
       'bottomMargin', 0.1, ...
       'topMargin', 0.05);
   
   
    ax = subplot('Position', subplotPosVectors(1,1).v);
    tmp = subplotPosVectors(1,2).v;
    dd = 0.02;
    tmp(2) = tmp(2) + dd;
    tmp(4) = tmp(4) - 2*dd;
    ax2 = subplot('Position', tmp);
    
    tmp = subplotPosVectors(1,3).v;
    tmp(2) = tmp(2) + dd;
    tmp(4) = tmp(4) - 2*dd;
    ax3 = subplot('Position', tmp);

    % Visualized best looking trials
    visualizedTrials = 1:nTrials;
    visualizedTrials = [3 5 6 7]; % 8 9 10 11 12 13 14 15 16];
    
    filtersNum = numel(videoData.filterBank);
    filterColors = brewermap(filtersNum, 'spectral');
    rowColPosition = zeros(filtersNum,2);
    filterRows = videoData.filterBank{1}.rowsNum;
    filterCols = videoData.filterBank{1}.colsNum;
    for k = 1:filtersNum
        rowColPosition(k,:) = videoData.filterBank{k}.rowColPosition; 
    end
    
    micronsToArcMin = 60/theMosaic.micronsPerDegree;
    emPathArcMin = videoData.stimResponseInstances.theMosaicEyeMovementsMicrons*micronsToArcMin;
    emRange = max(abs(emPathArcMin(:)))*[-1 1];
    deltaEM = emRange(2)/3;
    
    emTick = emRange(1):deltaEM:emRange(2);
    emTickLabel = sprintf('%2.1f\n', emTick);
            
    for ii = 1:numel(visualizedTrials)
        % Get trial data
        iTrial = visualizedTrials(ii);
        theDifferentialResponses = squeeze(videoData.differentialResponses(iTrial,:,:));
        theEMPath = squeeze(emPathArcMin(iTrial,:,:));
    
        if (iscell(videoData.filterBank))
            filterResponse = zeros(filtersNum, numel(tBinIndices));
            for k = 1:filtersNum
                % compute response of filter
                filterResponse(k,:) = computeFilterResponse(videoData.filterBank{k}, theDifferentialResponses, tBinIndices);
            end
        else
            filterResponse(1,:) = computeFilterResponse(videoData.filterBank, theDifferentialResponses, tBinIndices);
        end
        
        % Render a video frame for each time bin
        for t = 1:numel(tBinIndices)
            tt = tBinIndices(t);
            % The mosaic response
            instantaneousDiffResponse = theDifferentialResponses(:,tt);
            theMosaic.renderActivationMap(ax, instantaneousDiffResponse , ...
                    'visualizedConeAperture', 'geometricArea', ...
                    'mapType', 'modulated disks', ...
                    'signalRange', responseRange, ...
                    'colorMap', responseColorMap, ...
                    'showColorBar', ~true, ...
                    'labelColorBarTicks', ~true, ...
                    'outlineConesAlongHorizontalMeridian', ~true, ...
                    'showXLabel', ~true, ...
                    'showYLabel', ~true, ...
                    'showXTicks', true, ...
                    'showYTicks', true, ...
                    'tickInc', 0.1, ...
                    'backgroundColor', 0*[0.5 0.5 0.5]);
            xlabel(ax, '\it space(degs)');
            % Superimpose pooling filter
            hold(ax, 'on');
            zLevels = [0.15 0.15];
            if (iscell(videoData.filterBank))
                % Ensemble of spatial pooling filters
                filtersNum = numel(videoData.filterBank);
                
                for k = 1:filtersNum
                    % display filter
                    filterColor = squeeze(filterColors(k,:));
                    xAxisMeters = videoData.filterBank{k}.xaxisDegs * theMosaic.micronsPerDegree * 1e-6;
                    yAxisMeters = videoData.filterBank{k}.yaxisDegs * theMosaic.micronsPerDegree * 1e-6;
            
                    contour(ax, xAxisMeters, yAxisMeters, videoData.filterBank{k}.RFprofile, zLevels, ...
                    'Color', [0 0 0], 'LineWidth', 6, 'LineStyle', '-');
                
                    contour(ax, xAxisMeters, yAxisMeters, videoData.filterBank{k}.RFprofile, zLevels, ...
                    'Color', filterColor, 'LineWidth', 4, 'LineStyle', '-');
                
                    % display filter temporal responses
                    plotFilterResponse(ax2, t, videoData.nullResponseInstances.timeAxis*1000, squeeze(filterResponse(k,:)), filterColor);
                    if (k == filtersNum)
                        hold(ax2, 'off');
                    end
                end
                % Plot instantaneous filter activations
                ax4 = axes('Position',[0.473 0.645 0.25 0.25]);
                plotInstantaneousFilterEnsembleActivation(ax4, t, filterResponse, rowColPosition, filterRows, filterCols, filterColors);
                    
            else
                % Single spatial pooling filter
                xAxisMeters = videoData.filterBank.xaxisDegs * theMosaic.micronsPerDegree * 1e-6;
                yAxisMeters = videoData.filterBank.yaxisDegs * theMosaic.micronsPerDegree * 1e-6;
            
                contour(ax, xAxisMeters, yAxisMeters, videoData.filterBank.RFprofile, zLevels, ...
                    'Color', 'r', 'LineWidth', 2, 'LineStyle', '-');
                
                % display filter response
                plotFilterResponse(ax2, t, videoData.nullResponseInstances.timeAxis*1000, filterResponse, filterColor);
                hold(ax2, 'off');
            end
            
            % Superimpose stimulus outline
            %plot(ax,theEMPathMeters(tt, 1)+videoData.stimulusOutlineMetersX, -theEMPathMeters(tt, 2)+videoData.stimulusOutlineMetersY, 'k-', 'LineWidth', 5);
            %plot(ax,theEMPathMeters(tt, 1)+videoData.stimulusOutlineMetersX, -theEMPathMeters(tt, 2)+videoData.stimulusOutlineMetersY, 'c-', 'LineWidth', 3);
            hold(ax, 'off');
            % Set visible space
            set(ax, 'XLim', visualizedMosaicRange, 'YLim', visualizedMosaicRange);
            title(ax,sprintf('mosaic activation (t:%2.0fms, n=%2.0f)', videoData.nullResponseInstances.timeAxis(t)*1000, ii), 'FontWeight', 'normal');
            
            yMax = max(abs(filterResponse(:)));
            dY = yMax/5;
            set(ax2, 'XLim', [videoData.nullResponseInstances.timeAxis(1) videoData.nullResponseInstances.timeAxis(end)]*1000, ...
                     'YLim', yMax * [0 1], 'FontSize', 18, 'YColor', 'none', ...
                     'XTick', 0:50:500, 'YTick', 0:dY:yMax);
            box(ax2, 'off');
            grid(ax2, 'on'); axis(ax2, 'square');
            xlabel(ax2, '\it time (msec)');
            title(ax2, 'pooling mechanisms activations', 'FontWeight', 'normal');
            
            % the current empath  
            plot(ax3, [-100 100], [0 0 ], 'k-'); hold(ax3, 'on');
            plot(ax3, [0 0 ], [-100 100], 'k-'); 
            plot(ax3,theEMPath(1:t,1), -theEMPath(1:t,2), 'k-', 'LineWidth', 4.0);
            plot(ax3,theEMPath(1:t,1), -theEMPath(1:t,2), 'r-', 'LineWidth', 2.0);
            plot(ax3, theEMPath(t,1), -theEMPath(t,2), 'rs', 'LineWidth', 1.5, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', 12);
            hold(ax3, 'off');
            set(ax3, 'XLim', emRange, 'YLim', emRange, 'XTick', emTick, 'YTick', emTick, ...
                'XTickLabel', emTickLabel, 'YTickLabel', emTickLabel, 'FontSize', 18);
            grid(ax3, 'on'); box(ax3, 'off'); axis(ax3, 'square');
            xlabel(ax3, '\it position (arc min)');
            ylabel(ax2, '\it position (arc min)');
            title(ax3, 'eye movement path', 'FontWeight', 'normal');
            drawnow;
            
            % Add video frame
            videoOBJ.writeVideo(getframe(hFig));
        end
    end
    
    % Close video stream
    videoOBJ.close();
    fprintf('File saved in %s\n', videoFileName);

end

function plotInstantaneousFilterEnsembleActivation(ax, t, filterResponse, rowColPosition, filterRows, filterCols, filterColor)
    outline.x = [-1 -1 1 1 -1]*0.5;
    outline.y = [-1 1 1 -1 -1]*0.5;
    maxFilterResponse = max(filterResponse(:));
    minFilterResponse = min(filterResponse(:));
    for k = 1:size(rowColPosition,1)
        gain = (filterResponse(k,t)-minFilterResponse)/(maxFilterResponse-minFilterResponse)*0.9 + 0.1;
        xx = outline.x*gain + rowColPosition(k,2);
        yy = outline.y*gain + rowColPosition(k,1);
        patch(ax,xx,yy, squeeze(filterColor(k,:)));
        hold(ax, 'on')
    end
    hold(ax, 'off')
    axis(ax, 'square'); box(ax, 'on'); grid(ax, 'on');
    backgroundColor = [1 1 1]; % [0.3 0.3 0.3]
    set(ax, 'XLim', (max(rowColPosition(:))+0.5)*[-1 1], 'YLim', (max(rowColPosition(:))+0.5)*[-1 1], 'Color', backgroundColor, ...
        'XTickLabel', {}, 'YTickLabel', {});
    
end

function plotFilterResponse(ax, t, timeAxis, filterResponse, filterColor)
    %rowColPosition
    plot(ax, timeAxis(1:t), filterResponse(1:t), 'k-', 'LineWidth', 4.0);
    hold(ax, 'on');
    plot(ax, timeAxis(1:t), filterResponse(1:t), 'k-', 'LineWidth', 2.0, 'Color', filterColor);          
end


function filterResponse = computeFilterResponse(filterBank, mosaicResponse, tBinIndices)
    mosaicResponse = (mosaicResponse(:,tBinIndices))';
    spatialDimension = 2;
    cosFilterLinearActivation = squeeze(sum(bsxfun(@times, mosaicResponse, filterBank.cosPhasePoolingWeights), spatialDimension));
    sinFilterLinearActivation = squeeze(sum(bsxfun(@times, mosaicResponse, filterBank.sinPhasePoolingWeights), spatialDimension));
    if strcmp(filterBank.activationFunction, 'linear')
            filterResponse = cosFilterLinearActivation;
        elseif strcmp(filterBank.activationFunction, 'energy')
            filterResponse = sqrt(cosFilterLinearActivation.^2 + sinFilterLinearActivation.^2);
        elseif (strcmp(filterBank.activationFunction,'fullWaveRectifier'))
            filterResponse  = abs(cosFilterLinearActivation) + abs(sinFilterLinearActivation);
        else
            error('Activation function (''%s''), must be either ''energy'' OR ''fullWaveRectifier''\n', filterBank.activationFunction);
    end
    
    filterResponse = filterResponse';
end

function videoData = extractVideoData(stimData, noStimData, filterBank, theMosaic, spatialParams, thresholdParams, signalSource)
    stimData
    stimData.responseInstanceArray
    videoData.stimResponseInstances = stimData.responseInstanceArray;
    
    noStimData
    noStimData.responseInstanceArray
    videoData.nullResponseInstances = noStimData.responseInstanceArray;
    
    videoData.filterBank = filterBank;
    
%    filterBank
    
%      struct with fields:
% 
%     cosPhasePoolingProfile: [512×512 double]
%     sinPhasePoolingProfile: [512×512 double]
%                  RFprofile: [512×512 double]
%                  xaxisDegs: [1×512 double]
%                  yaxisDegs: [1×512 double]
%     cosPhasePoolingWeights: [1×3706 double]
%     sinPhasePoolingWeights: [1×3706 double]
%     envelopePoolingWeights: [1×3706 double]
%                       type: 'V1QuadraturePair'
%         activationFunction: 'energy'
%          temporalPCAcoeffs: Inf
    
%     filterBank{1}
%     
%     
%      struct with fields:
% 
%                cyclesPerRF: 4
%             bandwidthIndex: 1
%                    rowsNum: 5
%                    colsNum: 5
%           orientationIndex: 1
%              orientationRF: 0
%                  ft2Dindex: 1
%             rowColPosition: [-2 -2]
%            spatialPosition: [-0.0263 -0.0263]
%     cosPhasePoolingProfile: [512×512 double]
%     sinPhasePoolingProfile: [512×512 double]
%                  RFprofile: [512×512 double]
%                  xaxisDegs: [1×512 double]
%                  yaxisDegs: [1×512 double]
%     cosPhasePoolingWeights: [1×3706 double]
%     sinPhasePoolingWeights: [1×3706 double]
%     envelopePoolingWeights: [1×3706 double]
%                       type: 'V1QuadraturePair'
%         activationFunction: 'energy'
%          temporalPCAcoeffs: Inf
         
    % The differential response instances
    if (strcmp(signalSource, 'isomerizations'))
        % subtract the mean isomerization response of each cone (estimated from the null stimulus)
        meanNoiseFreeIsomerizations = mean(noStimData.noiseFreeIsomerizations,2);
        meanNoiseFreeIsomerizations = reshape(meanNoiseFreeIsomerizations, [1 size(meanNoiseFreeIsomerizations,1) size(meanNoiseFreeIsomerizations,2)]);
        videoData.differentialResponses = bsxfun(@minus, videoData.stimResponseInstances.theMosaicIsomerizations,meanNoiseFreeIsomerizations);
    elseif (strcmp(signalSource, 'photocurrents'))
        % subtract the mean photocurrent response of each cone (estimated from the null stimulus)
        meanNoiseFreePhotocurrents = mean(noStimData.noiseFreePhotocurrents,2);
        meanNoiseFreePhotocurrents = reshape(meanNoiseFreePhotocurrents, [1 size(meanNoiseFreePhotocurrents,1) size(meanNoiseFreePhotocurrents,2)]);
        videoData.differentialResponses = bsxfun(@minus, videoData.stimResponseInstances.theMosaicPhotocurrents,meanNoiseFreePhotocurrents);
    else
        error('signalSource must be either ''isomerizations'' or ''photocurrents''.');
    end
    
    % The stimulus outline
    stimSizeDegs = spatialParams.gaussianFWHMDegs; % spatialParams.fieldOfViewDegs or spatialParams.gaussianFWHMDegs
    videoData.stimulusOutlineMetersX = 0.5*[-1 -1 1 1 -1] * stimSizeDegs * theMosaic.micronsPerDegree * 1e-6;
    videoData.stimulusOutlineMetersY = 0.5*[-1 1 1 -1 -1] * stimSizeDegs * theMosaic.micronsPerDegree * 1e-6;
end

function [stimData, noStimData] = loadResponseData(colorModulationParams, constantParamsList)

    % Get out some data we'll want
    readProgram = 't_coneCurrentEyeMovementsResponseInstances';
    colorModulationParamsTemp = colorModulationParams;
    colorModulationParamsTemp.coneContrasts = [0 0 0]';
    colorModulationParamsTemp.contrast = 0;
    paramsList = constantParamsList;
    paramsList{numel(paramsList)+1} = colorModulationParamsTemp;
    clear 'colorModulationParamsTemp'
    
    rwObject = IBIOColorDetectReadWriteBasic;
    ancillaryData = rwObject.read('ancillaryData',paramsList,readProgram);
    noStimData = rwObject.read('responseInstances',paramsList,readProgram);
    testConeContrasts = ancillaryData.testConeContrasts;
    testContrasts = ancillaryData.testContrasts;
    
    % Determine visualized condition index (index of condition with max testContrast)
    maxTestContrast = max(testContrasts);
    visualizedConditionIndex = [];
    parforConditionStructs = responseGenerationParforConditionStructsGenerate(testConeContrasts,testContrasts);
    for kk = 1:length(parforConditionStructs)
        thisConditionStruct = parforConditionStructs{kk};
        colorModulationParamsTemp = colorModulationParams;
        colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
        colorModulationParamsTemp.contrast = thisConditionStruct.contrast;
        if (colorModulationParamsTemp.contrast  == maxTestContrast)
            visualizedConditionIndex = kk;
        end
        paramsList = constantParamsList;
        paramsList{numel(paramsList)+1} = colorModulationParamsTemp;
        thisConditionStruct.paramsList = paramsList;
        parforConditionStructs{kk} = thisConditionStruct;
    end
 
    % Return stimData for the visualizedConditionIndex
    theVisualizedConditionStruct = parforConditionStructs{visualizedConditionIndex};
    paramsList = theVisualizedConditionStruct.paramsList;
    stimData = rwObject.read('responseInstances',paramsList,readProgram);
end

function theMosaic = loadTheMosaic(topLevelDirParams, mosaicParams)
    coneParamsList = {topLevelDirParams, mosaicParams};
    theProgram = 't_coneCurrentEyeMovementsResponseInstances';
    rwObject = IBIOColorDetectReadWriteBasic;
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
end


