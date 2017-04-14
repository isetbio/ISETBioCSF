function hFig = visualizeResponsesInstancesAndNoiseFreeResponsesAsDensityPlots(timeAxis, noStimResponseInstances, stimResponseInstances, noStimNoiseFreeResponse, stimNoiseFreeResponse,  p, responseLevels, responseLevelsNum, plotType, yAxisLabel, figNo)

    allInstances = [...
        noStimResponseInstances(:);...
        stimResponseInstances(:)...
        ];
    
    %responseLevels = prctile(allInstances(:), [p 100-p]);
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
                   'heightMargin',   0.05, ...
                   'widthMargin',    0.04, ...
                   'leftMargin',     0.09, ...
                   'rightMargin',    0.01, ...
                   'bottomMargin',   0.07, ...
                   'topMargin',      0.04);
        else
            responseNums = 1;
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                   'rowsNum', 1, ...
                   'colsNum', 2, ...
                   'heightMargin',   0.05, ...
                   'widthMargin',    0.04, ...
                   'leftMargin',     0.09, ...
                   'rightMargin',    0.01, ...
                   'bottomMargin',   0.13, ...
                   'topMargin',      0.07);
        end
    else
       if (ndims(noStimResponseInstances) == 2)
            responseNums = size(noStimResponseInstances,1);
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                   'rowsNum', 3, ...
                   'colsNum', 2, ...
                   'heightMargin',   0.05, ...
                   'widthMargin',    0.04, ...
                   'leftMargin',     0.09, ...
                   'rightMargin',    0.01, ...
                   'bottomMargin',   0.07, ...
                   'topMargin',      0.04);
        else
            responseNums = 1;
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                   'rowsNum', 1, ...
                   'colsNum', 2, ...
                   'heightMargin',   0.05, ...
                   'widthMargin',    0.04, ...
                   'leftMargin',     0.09, ...
                   'rightMargin',    0.01, ...
                   'bottomMargin',   0.13, ...
                   'topMargin',      0.07);
        end 
    end
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 700 50+responseNums*300], 'Color', [1 1 1]);
            
    for responseSegment = 1:2
        for responseIndex = 1:responseNums
            if (responseNums == 1)
                if (responseSegment == 1)
                    responseInstances = noStimResponseInstances;
                    noiseFreeResponse = noStimNoiseFreeResponse;
                    titlePostfix = 'null';
                else
                    responseInstances = stimResponseInstances;
                    noiseFreeResponse = stimNoiseFreeResponse;
                    titlePostfix = 'test';
                end
            else
                if (responseSegment == 1)
                    responseInstances = squeeze(noStimResponseInstances(responseIndex,:,:));
                    noiseFreeResponse = squeeze(noStimNoiseFreeResponse(responseIndex,:));
                    titlePostfix = 'null';
                else
                    responseInstances = squeeze(stimResponseInstances(responseIndex,:,:));
                    noiseFreeResponse = squeeze(stimNoiseFreeResponse(responseIndex,:));
                    titlePostfix = 'test';
                end
            end
            % NO STIMULUS
            subplot('Position', subplotPosVectors(responseIndex,responseSegment).v);
            
            if (strcmp(plotType, 'density'))
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
            elseif (strcmp(plotType, 'line'))
                hold on;
                for instanceIndex = 1:size(responseInstances,1)
                    y = squeeze(responseInstances(instanceIndex, :));
                    xflip = [timeAxis(1 : end - 1) fliplr(timeAxis)];
                    yflip = [y(1 : end - 1) fliplr(y)];
                    patch(xflip, yflip, [0.1 0.1 0.1], 'EdgeAlpha', 0.1, 'FaceColor', 'none');
                end 
            end
            plot(timeAxis, noiseFreeResponse, 'w-', 'LineWidth', 6.0);
            plot(timeAxis, noiseFreeResponse, 'r-', 'LineWidth', 2.0);
            if (responseNums==3)
                x = timeAxisLimits(end)*0.75;
                y = double(responseLevels(end) - (responseLevels(end)-responseLevels(1))*0.05);
                if (responseIndex == 1)
                    text(x,y, 0.0, 'L-cone', 'FontSize', 12, 'FontWeight', 'bold');
                elseif (responseIndex == 2)
                    text(x,y, 0.0, 'M-cone', 'FontSize', 12, 'FontWeight', 'bold');
                elseif (responseIndex == 3)
                    text(x,y, 0.0, 'S-cone', 'FontSize', 12, 'FontWeight', 'bold');
                end
            end
            hold off;
            set(gca, 'CLim', [0 1]);
            set(gca, 'YLim', [responseLevels(1) responseLevels(end)], 'XLim', timeAxisLimits, 'FontSize', 14);
            axis 'xy';
            if (responseSegment == 1)
                ylabel(yAxisLabel, 'FontWeight', 'bold');
            else
                 set(gca, 'YTickLabel', {});
%                 posBefore = get(gca, 'Position');
%                 hcb = colorbar('eastoutside');
%                 colorTitleHandle = get(hcb,'Title');
%                 set(colorTitleHandle ,'String', sprintf(''));
%                 set(hcb, 'FontSize', 12);
%                 set(gca, 'Position', posBefore);
            end
            
            if (responseIndex == responseNums)
                xlabel('time (ms)', 'FontWeight', 'bold');
            end
            
            if (responseNums == 3)
                if (responseIndex == 1)
                    title(sprintf('%s', titlePostfix));
                end
            else
                title(titlePostfix);
            end
            grid on
            box on
        end % responseSegment
    end % responseIndex
    
    colormap(1-gray(1024));
    set(gcf,'renderer','opengl');
    drawnow;
end

