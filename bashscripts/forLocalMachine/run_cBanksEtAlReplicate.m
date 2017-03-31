function run_cBanksEtAlReplicate()

    % Original Banks et al mosaic params
    coneSpacingMicrons = 3.0;
    innerSegmentDiameter = 3.0;
    
    % Default isetbio mosaic params
    coneSpacingMicrons = 2.0;
    innerSegmentDiameter = 1.5797; % this corresponds to the 1.4 micron square pixel 
    
    nTrainingSamples =  128;
    cyclesPerDegree =  [10 20];
    luminances =  [3.4 34];
    nContrastsPerDirection =  12;
    lowContrast = 0.001;
    highContrast = 0.1;
        
    c_BanksEtAlReplicate(...
        'nTrainingSamples', nTrainingSamples, ...
        'cyclesPerDegree', cyclesPerDegree, ...
        'luminances', luminances, ...
        'nContrastsPerDirection', nContrastsPerDirection, ...
        'lowContrast', lowContrast, ...
        'highContrast', highContrast, ...
        'coneSpacingMicrons', coneSpacingMicrons, ...
        'innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(innerSegmentDiameter) ...
    );


end

