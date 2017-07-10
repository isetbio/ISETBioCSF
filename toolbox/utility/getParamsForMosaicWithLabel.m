function params = getParamsForMosaicWithLabel(mosaicName)
    switch mosaicName
        case 'Banks87'
            params.coneSpacingMicrons = 3.0;
            params.innerSegmentDiameter = 3.0;    % for a circular sensor
            params.conePacking = 'hexReg';
            params.LMSRatio = [0.67 0.33 0];
            params.mosaicRotationDegs = 30;
        case 'ISETbioHexRegNoScones'
            params.coneSpacingMicrons = 2.0;
            params.innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            params.conePacking = 'hexReg';
            params.LMSRatio = [0.67 0.33 0.0];
            params.mosaicRotationDegs = 0;
        case 'ISETbioHexRegLMS'
            params.coneSpacingMicrons = 2.0;
            params.innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            params.conePacking = 'hexReg';
            params.LMSRatio = [0.62 0.31 0.07];
            params.mosaicRotationDegs = 0;
        case 'ISETbioHexEccBasedNoScones'
            params.coneSpacingMicrons = 2.0;
            params.innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            params.conePacking = 'hex';
            params.LMSRatio = [0.67 0.33 0.0];
            params.mosaicRotationDegs = 0;
        case 'ISETbioHexEccBasedLMS'
            params.coneSpacingMicrons = 2.0;
            params.innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            params.conePacking = 'hex';
            params.LMSRatio = [0.62 0.31 0.07];
            params.mosaicRotationDegs = 0;
        otherwise
            error('Unknown mosaic: ''%s''.', mosaicName);
    end
end