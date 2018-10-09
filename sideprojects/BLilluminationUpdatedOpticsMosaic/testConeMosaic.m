function testConeMosaic
    
    resamplingFactor = 7;
    c = coneMosaicHex(resamplingFactor, ...
         'fovDegs', [0.4 0.4], ...
         'eccBasedConeDensity', true, ...
         'maxGridAdjustmentIterations', 200);
    % Visualize the cone mosaic and superimposed the achieved cone density map
    c.visualizeGrid('overlayConeDensityContour','measured');
  
    pause
    
    c = coneMosaic;
    fprintf('Default mosaic (rectangular): cone aperture area: %2.2f um^2\n', c.pigment.pdArea * 1e12);
    fprintf('Default mosaic (rectangular): cone aperture dimensions: %2.2f um x %2.2f um\n', c.pigment.pdWidth * 1e6, c.pigment.pdHeight * 1e6);
    fprintf('Default mosaic (rectangular): equivalent cone aperture diameter: %2.2f um\n', diameterForCircularApertureFromWidthForSquareAperture(c.pigment.pdWidth* 1e6));
    fprintf('Default mosaic (rectangular): cone separation: %2.2f um\n', c.pigment.width * 1e6);
    
    c = coneMosaic;
    wavelengthAxis = c.wave;
    quantalEfficiencies = c.qe;
    
    figure(1); clf;
    plot(wavelengthAxis, quantalEfficiencies(:,1), 'r-');
    hold on;
    plot(wavelengthAxis, quantalEfficiencies(:,2), 'g-');
    plot(wavelengthAxis, quantalEfficiencies(:,3), 'b-');
    ylabel('quantal efficiency');
    xlabel('wavelength');
    title('normal quantal efficiency');
    
    figure(2); clf;
    absorbances = c.pigment.absorbance;
    plot(wavelengthAxis, absorbances(:,1), 'r-');
    hold on;
    plot(wavelengthAxis, absorbances(:,2), 'g-');
    plot(wavelengthAxis, absorbances(:,3), 'b-');
    ylabel('absorbance');
    xlabel('wavelength');
    title('normal photopigment');
    
    % Make photopigment with anomalous absorbances
    absorbances(:,1) = 0.5*absorbances(:,1) + 0.5*absorbances(:,2);
    absorbances(:,2) = absorbances(:,2)*0;
    anomalousPhotoPigment = photoPigment('absorbance', absorbances);
    
    % Generate mosaic with anomalous photopigment
    c = coneMosaic('pigment', anomalousPhotoPigment);

    quantalEfficiencies = c.qe;
    figure(3); clf;
    plot(wavelengthAxis, quantalEfficiencies(:,1), 'r-');
    hold on;
    plot(wavelengthAxis, quantalEfficiencies(:,2), 'g-');
    plot(wavelengthAxis, quantalEfficiencies(:,3), 'b-');
    ylabel('quantal efficiency');
    xlabel('wavelength');
    title('anomalous quantal efficiency');
    
    absorbances = c.pigment.absorbance;
    figure(4); clf;
    plot(wavelengthAxis, absorbances(:,1), 'r-');
    hold on;
    plot(wavelengthAxis, absorbances(:,2), 'g-');
    plot(wavelengthAxis, absorbances(:,3), 'b-');
    ylabel('absorbance');
    xlabel('wavelength');
    title('anomalous photopigment');
    pause
    
    % Specify a mosaic at 0.5 degs (horizontal), 1.0 degs (vertical)
    eccentricityDegs = [0.5 1.0];
    eccentricityMicrons = eccentricityDegs * c.micronsPerDegree;
    eccentricityMeters = eccentricityMicrons * 1e-6;
    c = coneMosaic('center', eccentricityMeters);
    fprintf('Eccentric mosaic (rectangular): cone aperture area: %2.2f um^2\n', c.pigment.pdArea * 1e12);
    fprintf('Eccentric mosaic (rectangular): cone aperture dimensions: %2.2f um x %2.2f um\n', c.pigment.pdWidth * 1e6, c.pigment.pdHeight * 1e6);
    fprintf('Eccentric mosaic (rectangular): equivalent cone aperture diameter: %2.2f um\n', diameterForCircularApertureFromWidthForSquareAperture(c.pigment.pdWidth* 1e6));
    fprintf('Eccentric mosaic (rectangular): cone separation: %2.2f um\n', c.pigment.width * 1e6);
    
    
    fprintf('Hexagonal mosaic (they can only be placed at [0 0]\n')
    resamplingFactor = 7;
    c = coneMosaicHex(resamplingFactor, ...
        'fovDegs', [0.6 0.3]);
    fprintf('Default mosaic (rectangular): cone aperture area: %2.2f um^2\n', c.pigment.pdArea * 1e12);
    fprintf('Default mosaic (rectangular): cone aperture dimensions: %2.2f um x %2.2f um\n', c.pigment.pdWidth * 1e6, c.pigment.pdHeight * 1e6);
    fprintf('Default mosaic (rectangular): equivalent cone aperture diameter: %2.2f um\n', diameterForCircularApertureFromWidthForSquareAperture(c.pigment.pdWidth* 1e6));
    fprintf('Default mosaic (rectangular): cone separation: %2.2f um\n', c.pigment.width * 1e6);
    c.visualizeGrid();
    pause
    
    % But you can specify a desired aperture and separation:
    desiredConeSeparationInMicrons = 5.0;
    customApertureDiameterInMicrons = 3.0;
    c = coneMosaicHex(resamplingFactor, ...
        'fovDegs', [0.6 0.3], ...
        'customLambda', desiredConeSeparationInMicrons, ...
        'customInnerSegmentDiameter', customApertureDiameterInMicrons, ...
        'sConeFreeRadiusMicrons', 0);
    c.visualizeGrid();
    
    
    
end