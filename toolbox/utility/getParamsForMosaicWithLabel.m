function params = getParamsForMosaicWithLabel(mosaicName)

    params.sConeMinDistanceFactor = [];                     % no special treatment of S-cones
    params.sConeFreeRadiusMicrons = [];                     % no special treatment of S-cones
    params.latticeAdjustmentPositionalToleranceF = [];      
    params.latticeAdjustmentDelaunayToleranceF = [];
    params.marginF = [];
    params.eccBasedConeQuantalEfficiency = false;
    
    params.resamplingFactor = [];   % Empty indicates that c_BanksEtAlPhotocurrentAndEyeMovements
                                    % will choose a resamplingFactor based
                                    % on the mosaic size, smaller mosaics
                                    % will get higher resamplingFactor
    
    params.maxGridAdjustmentIterations = []; % Empty indicates that c_BanksEtAlPhotocurrentAndEyeMovements
                                    % will choose maxGridAdjustmentIterations based
                                    % on the mosaic size, smaller mosaics
                                    % will get higher rmaxGridAdjustmentIterations
                                    
    switch mosaicName
        case 'originalBanks'
            params.coneSpacingMicrons = 3.0;
            params.innerSegmentDiameter = 3.0;                      % used by Banks '87
            params.conePacking = 'hexReg';                          % regular hexagonal mosaic (no ecc variation)
            params.LMSRatio = [0.67 0.33 0];
            params.mosaicRotationDegs = 30;
        case 'originalBanksNoRotation'
            params.coneSpacingMicrons = 3.0;
            params.innerSegmentDiameter = 3.0;                      % used by Banks '87
            params.conePacking = 'hexReg';                          % regular hexagonal mosaic (no ecc variation)
            params.LMSRatio = [0.67 0.33 0];
            params.mosaicRotationDegs = 0;
        case 'ISETbioHexRegNoScones'
            params.coneSpacingMicrons = 2.0;
            params.innerSegmentDiameter = 1.5797;                   % for a circular sensor; this corresponds to ISETBio's default 1.4 micron square pixel (isetbio still does square apertures)
            params.conePacking = 'hexReg';                          % regular hexagonal mosaic (no ecc variation)
            params.LMSRatio = [0.67 0.33 0.0];
            params.mosaicRotationDegs = 0;
        case 'ISETbioHexRegLMS'
            params.coneSpacingMicrons = 2.0;
            params.innerSegmentDiameter = 1.5797;                   % for a circular sensor; this corresponds to ISETBio's default 1.4 micron square pixel (isetbio still does square apertures)
            params.conePacking = 'hexReg';                          % regular hexagonal mosaic (no ecc variation)
            params.LMSRatio = [0.62 0.31 0.07];
            params.mosaicRotationDegs = 0;
        case 'ISETbioHexEccBasedNoScones'
            params.coneSpacingMicrons = 2.0;
            params.innerSegmentDiameter = 1.5797;                   % for a circular sensor; this corresponds to ISETBio's default 1.4 micron square pixel (isetbio still does square apertures)
            params.conePacking = 'hex';                             % ecc-based cone density
            params.eccBasedConeQuantalEfficiency = false;           % no correction of cone efficiency with eccentricity
            params.LMSRatio = [0.67 0.33 0.0];
            params.mosaicRotationDegs = 0;
            params.sConeMinDistanceFactor = [];                     % no special treatment of S-cones
            params.sConeFreeRadiusMicrons = [];                     % no special treatment of S-cones
            params.latticeAdjustmentPositionalToleranceF = [];      
            params.latticeAdjustmentDelaunayToleranceF = [];
            params.marginF = [];                                    % marginF set according to FOV
        case 'ISETbioHexEccBasedLMS'
            params.coneSpacingMicrons = 2.0;
            params.innerSegmentDiameter = 1.5797;                   % for a circular sensor; this corresponds to ISETBio's default 1.4 micron square pixel (isetbio still does square apertures)
            params.conePacking = 'hex';                             % ecc-based cone density
            params.eccBasedConeQuantalEfficiency = false;           % no correction of cone efficiency with eccentricity
            params.LMSRatio = [0.62 0.31 0.07];
            params.mosaicRotationDegs = 0;
            params.sConeMinDistanceFactor = [];                     % no special treatment of S-cones
            params.sConeFreeRadiusMicrons = [];                     % no special treatment of S-cones
            params.latticeAdjustmentPositionalToleranceF = [];
            params.latticeAdjustmentDelaunayToleranceF = [];
            params.marginF = [];                                    % marginF set according to FOV
        case 'ISETbioHexEccBasedLMSrealistic'
            params.coneSpacingMicrons = 2.0;
            params.innerSegmentDiameter = 1.5797;                   % for a circular sensor; this corresponds to ISETBio's default 1.4 micron square pixel (isetbio still does square apertures)
            params.conePacking = 'hex';                             % ecc-based cone density
            params.eccBasedConeQuantalEfficiency = false;           % no correction of cone efficiency with eccentricity
            params.LMSRatio = [0.60 0.30 0.10];                     % More S-cones because S-cones in the S-cone free region will be eliminated - Also needed to generate a different folder
            params.mosaicRotationDegs = 0;
            params.sConeMinDistanceFactor = 3;                      % no special treatment of S-cones
            params.sConeFreeRadiusMicrons = 45;                     % no special treatment of S-cones
            params.latticeAdjustmentPositionalToleranceF = [];
            params.latticeAdjustmentDelaunayToleranceF = [];
            params.marginF = [];                                    % marginF set according to FOV
        case 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection'
            params.coneSpacingMicrons = 2.0;
            params.innerSegmentDiameter = 1.5797;                   % for a circular sensor; this corresponds to ISETBio's default 1.4 micron square pixel (isetbio still does square apertures)
            params.conePacking = 'hex';                             % ecc-based cone density
            params.eccBasedConeQuantalEfficiency = true;            % apply correction to cone efficiency with eccentricity
            params.LMSRatio = [0.60 0.30 0.10];                     % More S-cones because S-cones in the S-cone free region will be eliminated - Also needed to generate a different folder
            params.mosaicRotationDegs = 0;
            params.sConeMinDistanceFactor = 3;                      % no special treatment of S-cones
            params.sConeFreeRadiusMicrons = 45;                     % no special treatment of S-cones
            params.latticeAdjustmentPositionalToleranceF = [];
            params.latticeAdjustmentDelaunayToleranceF = [];
            params.marginF = [];                                    % marginF set according to FOV
        
        otherwise
            error('Unknown mosaic: ''%s''.', mosaicName);
    end
    
   
end