function debugCustomOI
    
    oiParams.pupilDiamMm = 3.5; 
    oiParams.umPerDegree = 300;
    oiParams.wavefrontSpatialSamples = 201;
    
    opticsModels = {...
        'AOoptics75mmPupil' ...
%         'WvfHumanSubject1' ...
%         'WvfHumanSubject2' ...
%         'WvfHumanSubject3' ...
%         'WvfHumanSubject4' ...
%         'WvfHumanSubject5' ...
%         'WvfHumanMeanOTFmagZeroOTFphase' ...
%         'WvfHumanMeanOTFmagMeanOTFphase' ...
        };
    
    PSFrange = 3;
    OTFrange = 100;

    for k = 1:numel(opticsModels)
        oiParams.opticsModel = opticsModels{k};
        [theCustomOI{k}, Zcoeffs] = oiWithCustomOptics(oiParams.opticsModel, oiParams.wavefrontSpatialSamples, oiParams.pupilDiamMm, oiParams.umPerDegree);
        visualizeThePSFs(theCustomOI{k}, oiParams.opticsModel, PSFrange, OTFrange, 'PSF', 100+k);
        visualizeThePSFs(theCustomOI{k}, oiParams.opticsModel, PSFrange, OTFrange, 'OTF', 200+k);
    end
    
end

function visualizeThePSFs(theOI, opticsModel, PSFrange, OTFrange, functionToVisualize, figNo)
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 4, ...
           'colsNum', 8, ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.02, ...
           'topMargin',      0.02);
       
    optics = oiGet(theOI,'optics');
    wavelengthsList = opticsGet(optics,'wave');
    fullOTF = opticsGet(optics, 'otf data');
    xSfGridCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
    ySfGridCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
    [xSfGridCyclesDegGrid,ySfGridCyclesDegGrid] = meshgrid(xSfGridCyclesDeg,ySfGridCyclesDeg);
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 2538 1289], 'Name', opticsModel);
    
    if (strcmp(functionToVisualize, 'PSF'))
        range = PSFrange;
    else
        range = OTFrange;
    end
    
    for waveIndex = 1:numel(wavelengthsList)
        theWaveOTF = fullOTF(:,:,waveIndex);
        [xGridMinutes,yGridMinutes,theWavePSF] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(theWaveOTF));
        row = floor((waveIndex-1)/size(subplotPosVectors,2))+1;
        col = mod(waveIndex-1, size(subplotPosVectors,2))+1;
        subplot('Position', subplotPosVectors(row,col).v);
        %contourLevels = 0:0.01:1.0;
        %contourf(xGridMinutes, yGridMinutes, theWavePSF/max(theWavePSF(:)), contourLevels, 'LineColor', 'none');
        if (strcmp(functionToVisualize, 'PSF'))
            h = pcolor(xGridMinutes, yGridMinutes, theWavePSF/max(theWavePSF(:)));
            xAx = xGridMinutes(1,:);
            yAx = yGridMinutes(:,1);
            centerPosition = floor(length(xAx)/2)+1;
            slice = theWavePSF(centerPosition,:)/max(theWavePSF(:));
        else
            otfMag = abs(fftshift(theWaveOTF));
            xAx = xSfGridCyclesDegGrid(1,:);
            yAx = ySfGridCyclesDegGrid(:,1);
            centerPosition = floor(length(xAx)/2)+1;
            h = pcolor(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid, otfMag);
            slice = otfMag(centerPosition,:);
        end
        set(h, 'EdgeColor', 'none');
    
        hold on;

        
        plot(xAx(centerPosition)*[1 1], [yAx(1) yAx(end)], 'k-');
        plot([xAx(1) xAx(end)], yAx(centerPosition)*[1 1], 'k-');
        plot(xAx, slice*range, 'rs-', 'LineWidth', 1.5);
        axis 'square'
        
        set(gca, 'XLim', range*[-1 1], 'YLim', range*[-1 1], 'CLim', [0 1]);
        if (row == size(subplotPosVectors,1)) && (col == 1)
            if (strcmp(functionToVisualize, 'PSF'))
                xlabel('min arc');
            else
                xlabel('c/deg');
            end
        else
            set(gca, 'XTickLabel', {}, 'YTickLabel', {});
        end
        title(sprintf('%d nm\nPSF volume: %1.5f, OTF(1,1): %1.5f',  wavelengthsList(waveIndex), sum(theWavePSF(:)), theWaveOTF(1,1)));
        drawnow;
        colormap(gray(1024)); 
    end
end
