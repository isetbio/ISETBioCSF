function [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs] = paramsForComparativeMosaicAnalysis(mosaicType)
    switch mosaicType
        case 'originalBanks'
            % 1. Original Banks et al mosaic params
            coneSpacingMicrons = 3.0;
            innerSegmentDiameter = 3.0;    % for a circular sensor
            conePacking = 'hexReg';
            LMSRatio = [0.67 0.33 0];
            mosaicRotationDegs = 30;
            
        case 'defaultIsetbio'
            % 2. Default isetbio mosaic params
            coneSpacingMicrons = 2.0;
            innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            conePacking = 'hexReg';
            LMSRatio = [0.60 0.30 0.10];
            mosaicRotationDegs = 0;
            
        case 'fullIsetbioNoScones'
            % 3. spatially-varying density isetbio mosaic params
            coneSpacingMicrons = 2.0;
            innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            conePacking = 'hex';
            LMSRatio = [0.67 0.33 0];
            mosaicRotationDegs = 0;
            
        case 'fullIsetbioWithScones'
            % 3. spatially-varying density isetbio mosaic params
            coneSpacingMicrons = 2.0;
            innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            conePacking = 'hex';
            LMSRatio = [0.60 0.30 0.10];
            mosaicRotationDegs = 0;
    end