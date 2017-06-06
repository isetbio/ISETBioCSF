function figGenerateOpticalFactorsBanksMosaic()

    loadParams();
    opticsModels = {'Geisler', 'WvfHuman'};
    
    for opticsModelIndex = 1:numel(opticsModels)
        [~,~,~, theOIs{opticsModelIndex}] = c_BanksEtAlPhotocurrentAndEyeMovements(...
            'opticsModel', opticsModels{opticsModelIndex}, ...
            'imagePixels', imagePixels, ...
            'cyclesPerDegree', cyclesPerDegreeExamined, ...
            'luminances', luminancesExamined, ...
            'nTrainingSamples', nTrainingSamples, ...
            'lowContrast', lowContrast, ...
            'highContrast', highContrast, ...
            'nContrastsPerDirection', nContrastsPerDirection, ...
            'ramPercentageEmployed', ramPercentageEmployed, ...
            'emPathType', emPathType, ...
            'centeredEMPaths', centeredEMPaths, ...
            'responseStabilizationMilliseconds', responseStabilizationMilliseconds, ...
            'responseExtinctionMilliseconds', responseExtinctionMilliseconds, ...
            'freezeNoise', freezeNoise, ...
            'integrationTime', integrationTimeMilliseconds/1000, ...
            'coneSpacingMicrons', coneSpacingMicrons, ...
            'innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(innerSegmentDiameter), ...
            'conePacking', conePacking, ...
            'LMSRatio', LMSRatio, ...
            'mosaicRotationDegs', mosaicRotationDegs, ...
            'computeMosaic', false, ...
            'computeOptics', false, ...
            'visualizeMosaic', false, ...
            'visualizeOptics', true, ...
            'computeResponses', false, ...
            'computePhotocurrentResponseInstances', false, ...
            'visualizeResponses', false, ...
            'visualizeSpatialScheme', false, ...
            'findPerformance', false, ...
            'visualizePerformance', false, ...
            'performanceSignal' , performanceSignal, ...
            'performanceClassifier', performanceClassifier ...
        );
    end
    
    visualizeOpticsModels(theOIs, opticsModels);
    
    function loadParams()
       imagePixels = 1024;
    
        % 'random'; 'frozen0';
        emPathType = 'frozen0'; %random'; %'random';     
        centeredEMPaths = false;

        % Use a subset of the trials. Specify [] to use all available trials
        nTrainingSamples = 1024;

        % Mosaic params
        coneSpacingMicrons = 3.0;
        innerSegmentDiameter = 3.0;    % for a circular sensor
        conePacking = 'hexReg';
        LMSRatio = [0.67 0.33 0];
        mosaicRotationDegs = 30;
        integrationTimeMilliseconds =  5.0;

        % response params
        responseStabilizationMilliseconds = 10;
        responseExtinctionMilliseconds = 50;

        % Conditions 
        lowContrast = 0.0001;
        highContrast = 0.3;
        nContrastsPerDirection =  18;
        luminancesExamined =  [34];

        % 'isomerizations', 'photocurrents'
        performanceSignal = 'isomerizations';

         %'mlpt'% 'svmV1FilterBank';
        performanceClassifier = 'mlpt';

        % Freeze noise for repeatable results
        freezeNoise = true;
        
        cyclesPerDegreeExamined = [2.5]; % [2.5 5 10 20 40 50];
        
        ramPercentageEmployed = 1.0;
    end

end

function visualizeOpticsModels(theOIlist, opticsModels)
    
    wavelengthsToCompute = 450:25:650;
    
    psfMax = 0;
    for waveIndex = 1:numel(wavelengthsToCompute) 
        theWavelength = wavelengthsToCompute(waveIndex);
        for k = 1:numel(theOIlist)
            theOI = theOIlist{k};
            [theOTF, otf_fx, otf_fy, thePSF, psf_x, psf_y] = getOtfPsfData(theOI, theWavelength);
            if (max(abs(thePSF(:))) > psfMax)
                psfMax = max(abs(thePSF(:)));
            end
        end
    end
    
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1]);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'colsNum', numel(wavelengthsToCompute), ...
           'rowsNum', 2+numel(theOIlist), ...
           'heightMargin',   0.08, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.04);
       
   
    sfRange = [0 80];
    psfRange = [0 1.0];
    psfTicks = [-1:0.2:1];
    
    modelColors = [1 0 0; 0 0 1; 0 0 0];
    
    for waveIndex = 1:numel(wavelengthsToCompute)  
        theWavelength = wavelengthsToCompute(waveIndex);
        for k = 1:numel(theOIlist)
            theOI = theOIlist{k};
            [theOTF, otf_fx, otf_fy, thePSF, psf_x, psf_y] = getOtfPsfData(theOI, theWavelength);
 
            selectedRow = floor(size(theOTF,1)/2) + 1;
            otfSlice = squeeze(theOTF(selectedRow,:));
            psfSlice = squeeze(thePSF(selectedRow,:));

            subplot('Position', subplotPosVectors(1,waveIndex).v);
            hold on;
            plot(otf_fx, otfSlice, '-', 'Color',  squeeze(modelColors(k,:)), 'LineWidth', 1.5);
            set(gca, 'XLim', sfRange, 'YLim', [0 1], 'FontSize', 12);
            if (waveIndex == 1)
                ylabel('OTF', 'FontWeight', 'bold');
            else
                set(gca, 'YTickLabels', {});
            end
            xlabel('c/deg', 'FontWeight', 'bold');
            axis 'square'
            title(sprintf('%d nm', theWavelength));
            
            subplot('Position', subplotPosVectors(2,waveIndex).v);
            hold on
            plot(psf_x, psfSlice, '-', 'Color',  squeeze(modelColors(k,:)), 'LineWidth', 1.5);
            set(gca, 'XTick', psfTicks, 'XLim', psfRange, 'YLim', [0 psfMax], 'FontSize', 12);
            if (waveIndex == 1)
                ylabel('PSF', 'FontWeight', 'bold');
            else
                set(gca, 'YTickLabels', {});
            end
            xlabel('arc min', 'FontWeight', 'bold');
            axis 'square'
            
            subplot('Position', subplotPosVectors(2+k,waveIndex).v);
            imagesc(psf_x, psf_y, thePSF/max(thePSF(:)));
            axis 'image'; axis 'xy';
            set(gca, 'XTick', psfTicks, 'YTick', psfTicks, 'XLim', psfRange, 'YLim', psfRange, 'CLim', [0 1], 'FontSize', 12);
            if (waveIndex == 1)
                ylabel(sprintf('%s PSF', opticsModels{k}), 'FontWeight', 'bold');
            end
            xlabel('arc min', 'FontWeight', 'bold');
            box on; grid on;
            colormap(gray(1024));
        end

        subplot('Position', subplotPosVectors(1,waveIndex).v);
        box on; grid on;
        legend(opticsModels);

        subplot('Position', subplotPosVectors(2,waveIndex).v);
        box on; grid on;
        legend(opticsModels);
    end
    
    
end

function [otf, otf_fx, otf_fy, psf, psf_x, psf_y] = getOtfPsfData(theOI, wavelengthsToCompute)
    theOptics = oiGet(theOI, 'optics');
    otf = abs(fftshift(opticsGet(theOptics,'otf data', wavelengthsToCompute)));
    otf_fx = opticsGet(theOptics, 'otf fx', 'cyclesperdeg');
    otf_fy = opticsGet(theOptics, 'otf fx', 'cyclesperdeg');

    psf = opticsGet(theOptics,'psf data', wavelengthsToCompute);
    psfSupport = opticsGet(theOptics,'psf support');
    psf_x = psfSupport{1}(1,:)*60;
    psf_y = psf_x;
end
