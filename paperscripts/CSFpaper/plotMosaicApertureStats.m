function plotMosaicApertureStats

    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 3;
    
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

    params.mosaicRotationDegs = 360
    params.cyclesPerDegreeExamined = [32 50 60];
    
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

    coneDensityLevels = 1e3 * [140 160 180 210 240 270];
    
    for sfIndex = 1:numel(params.cyclesPerDegreeExamined)
        theCurrentMosaic = theMosaicTypes.theMosaics{sfIndex};
        %theCurrentMosaic.displayInfo('plotApertureStats', true);
        
        theCurrentMosaic.visualizeGrid('generateNewFigure', true, ...
            'visualizedConeAperture', 'geometricArea', ...
            'overlayConeDensityContour', 'theoretical_and_measured', ...
            'coneDensityContourLevels', coneDensityLevels, ...
            'labelConeTypes', false, ...
            'overlayContourLabels', true, ...
            'BackgroundColor', [0.8 0.8 0.8]);
        
    end
end
