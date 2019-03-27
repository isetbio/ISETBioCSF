function testBleaching
    %load('/Volumes/DropBoxDisk/Dropbox/Dropbox (Aguirre-Brainard Lab)/IBIO_analysis/IBIOColorDetect/[c_BanksEtAlPhotocurrentAndEyeMovements]/M_hexPacking_coneSizeUm1.5797_coneSepUmNaN_VariedConeEff_VariedMP_rotationDegs360_eccDegs0.00_LMSdensities0.60_0.30_0.10_FOVdegs3.94x3.94_intTimeMs5_photonNoiseFrozen_osModelLinear_osTimeStepMs0.10_osNoiseFrozen_apBlur1_dark_0_0_0/matfiles/t_coneCurrentEyeMovementsResponseInstances/coneMosaic.mat')
    %theConeMosaic = theData;
    %theConeMosaic.displayInfo()
    
    innerSegmentAreaMicronsSquared(1) = 1.96;
    innerSegmentAreaMicronsSquared(2) = 9.2;
    innerSegmentAreaMicronsSquared(3) = 19.6785;
    
    % Intensity of steady background in photons * sec^{-1} * micron^{2}
    iHalf = 10^5.57;
    
    %                                   1%  50%   99%
    iBackgroundPhotonsPerSecPerCone = [317 14182 28136];
    
    % Compute background in photons * sec^{-1} * micron^{2}
    iBackground = iBackgroundPhotonsPerSecPerCone ./innerSegmentAreaMicronsSquared ;
    
    % Compute percent bleached
    % Intensity of steady background in photons * sec^{-1} * micron^{2},
    % from "Light adaptation and photopigment bleaching in cone
    % photoreceptors in situ in the retina of the turtle" by DA Burkhardt
    iHalf = 10^5.57;
    fractionPhotopigmentBleached = 1 - iHalf ./(iHalf + iBackground)
    
    % Reproduce Figure 4 of "Light adaptation and photopigment bleaching in cone
    % photoreceptors in situ in the retina of the turtle" by DA Burkhardt
    figure()
    iBackground = (0:1:10);
    percentPhotopigment  =  iHalf ./(iHalf + 10.^iBackground);
    plot(iBackground, percentPhotopigment, 'ks-');
end

