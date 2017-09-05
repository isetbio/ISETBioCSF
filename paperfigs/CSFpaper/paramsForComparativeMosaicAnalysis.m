function [coneSpacingMicrons, innerSegmentDiameter, conePacking, LMSRatio, mosaicRotationDegs, mosaicLegend] = paramsForComparativeMosaicAnalysis(mosaicType)
 
error('Call params = getParamsForMosaicWithLabel(mosaicName) instead');

    switch mosaicType
        case 'originalBanks'
            % 1. Original Banks et al mosaic params
            coneSpacingMicrons = 3.0;
            innerSegmentDiameter = 3.0;    % for a circular sensor
            conePacking = 'hexReg';
            LMSRatio = [0.67 0.33 0];
            mosaicRotationDegs = 30;
            mosaicLegend = 'ISETbio-Banks (LM,  hexReg@3.0um,  a=3.0um)';
            
        case 'originalBanksNoRotation'
            % 1. Original Banks et al mosaic params
            coneSpacingMicrons = 3.0;
            innerSegmentDiameter = 3.0;    % for a circular sensor
            conePacking = 'hexReg';
            LMSRatio = [0.67 0.33 0];
            mosaicRotationDegs = 0;
            mosaicLegend = 'ISETbio-Banks (LM,  hexReg@3.0um,  a=3.0um)';
            
        case 'defaultIsetbio'
            % 2. Default isetbio mosaic params
            coneSpacingMicrons = 2.0;
            innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            conePacking = 'hexReg';
            LMSRatio = [0.67 0.33 0.0];
            mosaicRotationDegs = 0;
            mosaicLegend = 'ISETbio-LMreg (LM,  hexReg@2.0um,  a=1.6um)';
            
        case 'fullIsetbioNoScones'
            % 3. spatially-varying density isetbio mosaic params
            coneSpacingMicrons = 2.0;
            innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            conePacking = 'hex';
            LMSRatio = [0.67 0.33 0];
            mosaicRotationDegs = 0;
            mosaicLegend = 'ISETbio-LM (LM,  hex@ecc-based, a=1.6um)';
            
        case 'fullIsetbioWithScones'
            % 3. spatially-varying density isetbio mosaic params
            coneSpacingMicrons = 2.0;
            innerSegmentDiameter = 1.5797; % for a circular sensor; this corresponds to the 1.4 micron square pixel 
            conePacking = 'hex';
            LMSRatio = [0.62 0.31 0.07];
            mosaicRotationDegs = 0;
            mosaicLegend = 'ISETbio-full  (LMS, hex@ecc-based, a=1.6um)';
    end