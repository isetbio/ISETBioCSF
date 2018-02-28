function visualizeAllCustomWVFOpticsModels()

    wavefrontSpatialSamples = 261*2+1;
    calcPupilDiameterMM = 3;
    umPerDegree = 300;
    psfRange = 6;
    showTranslation = false;
    
    opticsModels = availableCustomWvfOpticsModels();
    opticsModels = {opticsModels{2:end}};

    for oiIndex = 1:numel(opticsModels)
        
        [theCustomOI, Zcoeffs] = oiWithCustomOptics(opticsModels{oiIndex}, ...
            wavefrontSpatialSamples, calcPupilDiameterMM, umPerDegree, ...
            'showTranslation',showTranslation);
        
        optics = oiGet(theCustomOI, 'optics');
        wavelengths = opticsGet(optics, 'wave');
        wavelengthsToShow = 1:3:numel(wavelengths);
        
        if (oiIndex == 1)
            subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', numel(opticsModels), ...
                'colsNum', numel(wavelengthsToShow), ...
                'heightMargin', 0.03, ...
                'widthMargin', 0.02, ...
                'leftMargin', 0.03, ...
                'rightMargin', 0.001, ...
                'bottomMargin', 0.02, ...
                'topMargin', 0.02);
            hFig = figure(1); clf;
            set(hFig, 'Position', [10 10 1540 1345], 'Color', [1 1 1]);

        end

        xSfCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
        ySfCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
        [xSfGridCyclesDeg,ySfGridCyclesDeg] = meshgrid(xSfCyclesDeg,ySfCyclesDeg);
            
        for iW = 1:numel(wavelengthsToShow)
            targetWavelength = wavelengths(wavelengthsToShow(iW));
            
            waveOTF = opticsGet(optics,'otf data',targetWavelength);
            [xGridMinutes, yGridMinutes, wavePSF] = OtfToPsf(...
                xSfGridCyclesDeg,ySfGridCyclesDeg,fftshift(waveOTF), ...
                'warningInsteadOfErrorForNegativeValuedPSF', true ...
            );
            xMinutes = squeeze(xGridMinutes(1,:));
            yMinutes = squeeze(yGridMinutes(:,1));
            
            if (iW == 1)
                xx = find(abs(xMinutes) <= psfRange);
                yy = find(abs(yMinutes) <= psfRange);
                colormap(gray(1024));
            end
            
            wavePSF = wavePSF(yy,xx);
            xMinutes = xMinutes(xx);
            yMinutes = yMinutes(yy);
            
            pos = subplotPosVectors(oiIndex, iW).v;
            
            subplot('Position', pos);
            imagesc(xMinutes, yMinutes,wavePSF);
            hold on;
            plot([0 0], [xMinutes(1) xMinutes(end)], 'r-', 'LineWidth', 1.0);
            plot([xMinutes(1) xMinutes(end)], [0 0], 'r-', 'LineWidth', 1.0);
            hold off;
            axis 'equal';
            set(gca, 'XLim', psfRange*1.05*[-1 1], 'YLim', psfRange*1.05*[-1 1], 'FontSize', 14);
            set(gca, 'XTick', -10:2:10, 'YTick', -10:2:10);
            grid on; box on;
            
            if (iW > 1)
                set(gca, 'YTickLabel', {});
                title(sprintf('\n%2.0fnm', targetWavelength), 'FontSize', 10);
            else
                title(sprintf('%s\n%2.0fnm', strrep(opticsModels{oiIndex},'Thibos','') , targetWavelength), 'FontSize', 10);
                ylabel('arc min');
            end
            
            if (oiIndex < numel(opticsModels))
                set(gca, 'XTickLabel', {});
            else
                xlabel('arc min');
            end
            
            drawnow;
        end
    end
    NicePlot.exportFigToPNG('customOpticsModels', hFig, 300);
end