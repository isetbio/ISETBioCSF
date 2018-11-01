function plotMosaicApertureStats

    % How to split the computation
    % 0 (All mosaics), 1; (Largest mosaic), 2 (Second largest), 3 (all 2 largest)
    computationInstance = 3;
    
    opticsName = 'ThibosBestPSFSubject3MMPupil';
     
    %mosaicName = 'ISETbioHexEccBasedLMSrealisticEfficiencyCorrection';
    mosaicName = 'ISETbioHexEccBasedLMSrealistic';
    
    params = getCSFpaperDefaultParams(mosaicName, computationInstance);
    
    
    params.LMSRatio = [0.67 0.33 0];
    %params.LMSRatio = [0.60 0.30 0.10];
    
    
    params.opticsModel = opticsName;
       
    params.coneContrastDirection = 'L+M+S';

    % Response duration params
    params.frameRate = 10; %(1 frames)
    params.responseStabilizationMilliseconds = 40;
    params.responseExtinctionMilliseconds = 40;
    
    
    
    % Eye movement params
    params.emPathType = 'frozen0';
    params.centeredEMpaths = ~true;

    if (strcmp(mosaicName,'originalBanks'))
        params.mosaicRotationDegs = 30;
    else
        params.mosaicRotationDegs = 360;
    end
    params.cyclesPerDegreeExamined = [2 4 8 16 32 50 60];
    
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
        
        [innerSegmentCoverage(sfIndex), geometricCoverage(sfIndex)] = theCurrentMosaic.retinalCoverage();
        
%         theCurrentMosaic.visualizeGrid('generateNewFigure', true, ...
%             'visualizedConeAperture', 'geometricArea', ...
%             'overlayConeDensityContour', 'theoretical_and_measured', ...
%             'coneDensityContourLevels', coneDensityLevels, ...
%             'labelConeTypes', false, ...
%             'overlayContourLabels', true, ...
%             'BackgroundColor', [0.8 0.8 0.8]);
        
    end
    
    figure(1)
    plot(params.cyclesPerDegreeExamined, innerSegmentCoverage, 'rs');
end
