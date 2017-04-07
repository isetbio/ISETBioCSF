function hFig = visualizeResponsesInstancesAndNoiseFreeResponsesAsDensityPlots(timeAxis, noStimResponseInstances, stimResponseInstances, noStimNoiseFreeResponse, stimNoiseFreeResponse,  p, responseLevelsNum, yAxisLabel)

    allInstances = [...
        noStimResponseInstances(:);...
        stimResponseInstances(:)...
        ];
    
    responseLevels = prctile(allInstances(:), [p 100-p]);
    responseLevels = linspace(responseLevels(1), responseLevels(2), responseLevelsNum);
    timeAxisLimits = [timeAxis(1) timeAxis(end)];
    if (timeAxis(1) == timeAxis(end))
        timeAxisLimits = timeAxis(1) + [-0.5 0.5];
    end
    
    if (numel(timeAxis) > 1)
        if (ndims(noStimResponseInstances) == 3)
            responseNums = size(noStimResponseInstances,1);
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                   'rowsNum', 3, ...
                   'colsNum', 2, ...
                   'heightMargin',   0.09, ...
                   'widthMargin',    0.05, ...
                   'leftMargin',     0.07, ...
                   'rightMargin',    0.07, ...
                   'bottomMargin',   0.07, ...
                   'topMargin',      0.04);
        else
            responseNums = 1;
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                   'rowsNum', 1, ...
                   'colsNum', 2, ...
                   'heightMargin',   0.09, ...
                   'widthMargin',    0.05, ...
                   'leftMargin',     0.07, ...
                   'rightMargin',    0.07, ...
                   'bottomMargin',   0.13, ...
                   'topMargin',      0.07);
        end
    else
       if (ndims(noStimResponseInstances) == 2)
            responseNums = size(noStimResponseInstances,1);
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                   'rowsNum', 3, ...
                   'colsNum', 2, ...
                   'heightMargin',   0.09, ...
                   'widthMargin',    0.05, ...
                   'leftMargin',     0.07, ...
                   'rightMargin',    0.07, ...
                   'bottomMargin',   0.07, ...
                   'topMargin',      0.04);
        else
            responseNums = 1;
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                   'rowsNum', 1, ...
                   'colsNum', 2, ...
                   'heightMargin',   0.09, ...
                   'widthMargin',    0.05, ...
                   'leftMargin',     0.07, ...
                   'rightMargin',    0.07, ...
                   'bottomMargin',   0.13, ...
                   'topMargin',      0.07);
        end 
    end
    
    hFig = figure(5001); clf;
    set(hFig, 'Position', [10 10 900 50+responseNums*300], 'Color', [1 1 1]);
            
    for responseSegment = 1:2
        for responseIndex = 1:responseNums
            if (responseNums == 1)
                if (responseSegment == 1)
                    responseInstances = noStimResponseInstances;
                    noiseFreeResponse = noStimNoiseFreeResponse;
                    titlePostfix = 'no stimulus';
                else
                    responseInstances = stimResponseInstances;
                    noiseFreeResponse = stimNoiseFreeResponse;
                    titlePostfix = 'test stimulus';
                end
            else
                if (responseSegment == 1)
                    responseInstances = squeeze(noStimResponseInstances(responseIndex,:,:));
                    noiseFreeResponse = squeeze(noStimNoiseFreeResponse(responseIndex,:));
                    titlePostfix = 'no stimulus';
                else
                    responseInstances = squeeze(stimResponseInstances(responseIndex,:,:));
                    noiseFreeResponse = squeeze(stimNoiseFreeResponse(responseIndex,:));
                    titlePostfix = 'test stimulus';
                end
            end
            % NO STIMULUS
            subplot('Position', subplotPosVectors(responseIndex,responseSegment).v);
            responseDistribution = zeros(numel(responseLevels)-1, numel(timeAxis));

            if (numel(timeAxis) > 1)
                for tBin = 1:numel(timeAxis)
                    N = histcounts(squeeze(responseInstances(:,tBin)),responseLevels);
                    responseDistribution(:,tBin) = N';
                end
            else
                N = histcounts(responseInstances(:),responseLevels);
                responseDistribution(:,1) = N';
            end
            imagesc(timeAxis, responseLevels, responseDistribution/max(responseDistribution(:)));
            hold on;
            plot(timeAxis, noiseFreeResponse, 'w-', 'LineWidth', 7.0);
            plot(timeAxis, noiseFreeResponse, 'r-', 'LineWidth', 2.0);
            hold off;
            set(gca, 'CLim', [0 1]);
            set(gca, 'YLim', [responseLevels(1) responseLevels(end)], 'XLim', timeAxisLimits, 'FontSize', 14);
            axis 'xy';
            if (responseSegment == 1)
                ylabel(yAxisLabel, 'FontWeight', 'bold');
            else
                set(gca, 'YTickLabel', {});
                posBefore = get(gca, 'Position');
                hcb = colorbar('eastoutside');
                colorTitleHandle = get(hcb,'Title');
                set(colorTitleHandle ,'String', sprintf(''));
                set(hcb, 'FontSize', 12);
                set(gca, 'Position', posBefore);
            end
            
            if (responseIndex == responseNums)
                xlabel('time (ms)', 'FontWeight', 'bold');
            end
            
            if (responseNums == 3)
                if (responseIndex == 1)
                    title(sprintf('%s (max responding L-cone)', titlePostfix));
                elseif (responseIndex == 2)
                    title(sprintf('%s (max responding M-cone)', titlePostfix));
                elseif (responseIndex == 3)
                    title(sprintf('%s (max responding S-cone)', titlePostfix));
                end
            else
                title(titlePostfix);
            end
            grid on
            box on
        end % responseSegment
    end % responseIndex
    
    colormap(1-gray(1024));
    drawnow;
end
