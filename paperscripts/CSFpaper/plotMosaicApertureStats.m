function plotMosaicApertureStats

    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 0;
    
    opticsName = 'ThibosBestPSFSubject3MMPupil';
     
    mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection';
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    params.opticsModel = opticsName;
       
    params.coneContrastDirection = 'L+M+S';

    % Response duration params
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;

    % Eye movement params
    params.emPathType = 'frozen0';
    params.centeredEMpaths = ~true;

    % Simulation steps to perform
    params.computeMosaic = ~true; 
    params.visualizeMosaic = ~true;

    params.computeResponses = ~true;
    params.computePhotocurrentResponseInstances = ~true;
    params.visualizeMosaic = true;
    params.visualizeResponses = ~true;
    params.visualizeSpatialScheme = ~true;
    params.visualizeOIsequence = ~true;
    params.visualizeOptics = ~true;
    params.visualizeMosaicWithFirstEMpath = ~true;
    params.visualizeSpatialPoolingScheme = ~true;
    params.visualizeStimulusAndOpticalImage = ~true;
    params.visualizeDisplay = ~true;

    params.visualizeKernelTransformedSignals = ~true;
    params.findPerformance = ~true;
    params.visualizePerformance = ~true;
    params.deleteResponseInstances = ~true;

    [theMosaicTypes, thePsychometricFunctions, theFigData] = ...
        run_BanksPhotocurrentEyeMovementConditions(params);

    for sfIndex = 4:4 % 1:numel(params.cyclesPerDegreeExamined)
        theCurrentMosaic = theMosaicTypes.theMosaics{sfIndex};
        theCurrentMosaic.displayInfo('plotApertureStats', true);
        theCurrentMosaic.displayInfo();
        idx = find(theCurrentMosaic.pattern ==3);
        theCurrentMosaic.pattern(idx) = 2;
        theCurrentMosaic.reassignConeIdentities();
        theCurrentMosaic.visualizeGrid('visualizedConeAperture', 'geometricArea', 'BackgroundColor', [0 0 0]);
            
    end
end
