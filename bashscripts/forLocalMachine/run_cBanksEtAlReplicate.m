function run_cBanksEtAlReplicate()

    mosaicType = 'originalBanks';
    
    switch mosaicType
        case 'originalBanks'
            % 1. Original Banks et al mosaic params
            coneSpacingMicrons = 3.0;
            innerSegmentDiameter = 3.0;    % for a circular sensor
            conePacking = 'hexReg';
            
        case 'defaultIsetbio'
            % 2. Default isetbio mosaic params
            coneSpacingMicrons = 2.0;
            innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            conePacking = 'hexReg';
    end
    
    freezeNoise = ~true;
    nTrainingSamples = 256;
    cyclesPerDegreeExamined =  [40];
    luminancesExamined =  [34];
    nContrastsPerDirection =  12;
    lowContrast = 0.001;
    highContrast = 0.3;
        
    c_BanksEtAlReplicate(...
        'freezeNoise', freezeNoise, ...
        'mosaicRotationDegs', 30, ...
        'nTrainingSamples', nTrainingSamples, ...
        'cyclesPerDegree', cyclesPerDegreeExamined, ...
        'luminances', luminancesExamined, ...
        'nContrastsPerDirection', nContrastsPerDirection, ...
        'lowContrast', lowContrast, ...
        'highContrast', highContrast, ...
        'coneSpacingMicrons', coneSpacingMicrons, ...
        'innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(innerSegmentDiameter) ...
    );


end

