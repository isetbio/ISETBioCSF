function removePhotocurrentDataFromResponseInstancesFile()

    conditionDir = 'CM_L71_M71_S0_con0.05000_wls_380_4_780';
    updateFile(conditionDir);
    
end

function updateFile(conditionDir)
    dropboxDir = '/home/isetbio/Data';
    experimentDir = '[c_BanksEtAlPhotocurrentAndEyeMovements]';
    mosaicDir = 'M_hexRegPacking_coneSizeUm1.4000_coneSepUm2.0000_rotationDegs0_eccentricityDegs0.00_LMSdensities0.60_0.30_0.00_FOVdegs3.15x3.15_intTimeMs5_photonNoiseRandom_osModelLinear_osTimeStepMs0.10_osNoiseRandom_apBlur0_dark_0_0_0';
    opticsDir = 'O_blur1_lens1_pupilDiam2.0_WvfHuman';
    spatialDir = 'SS_SFcpd2.50_phaseDegs0.00_FOVdegs3.15_HalfCos_FWHMDegs1.50';
    backgroundDir = 'BG_lum50.00_lumFactor0.680_leakLum1.0';
    temporalDir = 'ST_rampTauMilliSecsNaN_stimDurMilliSecs100_analyzedResponseDurationMilliSecs280_withOffsetMilliSecs0_eyeMovementPathFrozen0';
    colorDir = 'LM_angles45_to_90_N1_contrasts0.00050_to_0.800_N18_trialsN1024'; 
    dataFile = fullfile(dropboxDir, experimentDir, mosaicDir, opticsDir, spatialDir, temporalDir, backgroundDir, colorDir, conditionDir, 'matfiles/t_coneCurrentEyeMovementsResponseInstances/responseInstances.mat');
    
    fprintf('Loading data for condition ''%s'' ...\n', conditionDir);
    s = whos('-file', dataFile);
    if (~isfield(s, 'name')) || ((isfield(s, 'name')) && (~strcmp(s.name, 'theData')))
        fprintf('''theData'' variable was not found in the imported file. Doing nothing.\n');
        return;
    end
    load(dataFile, 'theData');
    
    fprintf('Emptying photocurrents\n');
    theData.responseInstanceArray.theMosaicPhotocurrents = [];
    
    fprintf('Updating %s ...\n', dataFile);
    save(dataFile, 'theData', '-v7.3');
    
    fprintf('Done!\n');
end

