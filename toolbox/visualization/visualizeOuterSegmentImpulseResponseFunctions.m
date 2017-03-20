function visualizeOuterSegmentImpulseResponseFunctions(timeAxis, noStimImpulseResponses, stimImpulseResponses, paramsList)

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                   'rowsNum', 3, ...
                   'colsNum', 1, ...
                   'heightMargin',   0.05, ...
                   'widthMargin',    0.01, ...
                   'leftMargin',     0.12, ...
                   'rightMargin',    0.01, ...
                   'bottomMargin',   0.05, ...
                   'topMargin',      0.02);
    
    hFig = figure(222); clf;
    set(hFig, 'Position', [10 10 300 970], 'Color', [1 1 1]);
    
    for coneIndex = 1:3
        subplot('Position', subplotPosVectors(coneIndex,1).v);
        plot(timeAxis, squeeze(noStimImpulseResponses(:,coneIndex)), 'k-', 'LineWidth', 1.5);
        set(gca, 'FontSize', 12, 'XLim', [timeAxis(1) timeAxis(end)], 'YLim', [min(squeeze(noStimImpulseResponses(:,coneIndex))) max(squeeze(noStimImpulseResponses(:,coneIndex)))]);
        hold on;
        for kk = 1:size(stimImpulseResponses,1)
            plot(timeAxis, squeeze(stimImpulseResponses(kk,:,coneIndex)), 'k-', 'LineWidth', 1.5);
            set(gca, 'FontSize', 12, 'XLim', [timeAxis(1) timeAxis(end)], 'YLim', [min(squeeze(noStimImpulseResponses(:,coneIndex))) max(squeeze(noStimImpulseResponses(:,coneIndex)))]);
        end
        hold off;
        box off; grid on;
        if (coneIndex == 1)
            title('L-cone photocurrent impulse response');
        elseif (coneIndex == 2)
            title('M-cone photocurrent impulse response');
        elseif (coneIndex == 3)
            title('S-cone photocurrent impulse response');
        end
        if (coneIndex < 3)
            set(gca, 'XTickLabel', {});
        else
            xlabel('time (seconds)');
        end
    end
    drawnow;
    
    % Export figure  
    rwObject = IBIOColorDetectReadWriteBasic;
    theProgram = mfilename;      
    data = 0;
    rwObject.write('osImpulseResponseFunctions', data, paramsList, theProgram, ...
        'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
    
end

