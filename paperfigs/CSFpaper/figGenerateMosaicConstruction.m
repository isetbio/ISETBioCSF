function figGenerateMosaicConstruction()

%% Initialize
ieInit; clear; close all;

% Set random seed to obtain replicable results
rng(1235);

makeNew = true;
if (makeNew)
    % Set coneMosaicHex - specific params
    params.resamplingFactor = 9;                            % 9 is good;how fine to sample the hex mosaic positions with an underlying rect mosaic
    params.fovDegs = 0.5*[1 1];                             % FOV in degrees ([width height], default: 0.25x0.25
    params.eccBasedConeDensity = true;                      % if true, generate a mosaic where cone spacing varies with eccentricity (Curcio model)
    params.customLambda = [];                               % cone spacing in microns (only used with regular hex mosaics)
    params.rotationDegs = 0;                                % rotation of the mosaic, in degrees 0 = , 30 = (makes sense with regular hex mosaics)
    params.customInnerSegmentDiameter = [];                 % inner segment diameter, in microns (empty for isetbio default) 
    params.sConeMinDistanceFactor = 3.0;                    % min distance between neighboring S-cones = f * local cone separation - to make the S-cone lattice semi-regular
    params.sConeFreeRadiusMicrons = 45;                     % radius of S-cone free retina, in microns
    params.latticeAdjustmentPositionalToleranceF = 0.01;  % determines cone delta movement tolerance for terminating iterative adjustment - by default this is 0.01 (here setting it lower for faster, but less acurate mosaic generation)
    params.latticeAdjustmentDelaunayToleranceF = 0.001;   % determines position tolerance for triggering another Delaunay triangularization - by default this is 0.001 (here setting it lower for faster, but less acurate mosaic generation)
    params.marginF = 1.2;
    saveLatticeAdjustmentProgression = true;                % set to true, only if interested to see how the mosaic lattice is iteratively adjusted when eccBasedConeDensity is true               

    % Generate the mosaic
    theHexMosaic = coneMosaicHex(params.resamplingFactor, ...
        'fovDegs', params.fovDegs, ...
        'eccBasedConeDensity', params.eccBasedConeDensity, ...
        'rotationDegs', params.rotationDegs, ...
        'sConeMinDistanceFactor', params.sConeMinDistanceFactor, ...    
        'sConeFreeRadiusMicrons', params.sConeFreeRadiusMicrons, ...    
        'customLambda', params.customLambda, ...
        'customInnerSegmentDiameter', params.customInnerSegmentDiameter, ...
        'latticeAdjustmentPositionalToleranceF', params.latticeAdjustmentPositionalToleranceF, ...              
        'latticeAdjustmentDelaunayToleranceF', params.latticeAdjustmentDelaunayToleranceF, ...     
        'marginF', params.marginF, ...
        'saveLatticeAdjustmentProgression', saveLatticeAdjustmentProgression ...  
    );
    % Save the mosaic
    save('theHexMosaic.mat', 'theHexMosaic', '-v7.3');
else
    fprintf('\nLoading mosaic ... ');
    % Load the mosaic
    load('theHexMosaic1.5.mat');
    fprintf('Done\n');
end

% Display mosaic info(
theHexMosaic.displayInfo();

% Show how the lattice of an ecc-based cone density hex mosaic is iteratively adjusted
contourLevels = 1e3 * [75:25:250];
hFig = theHexMosaic.plotMosaicProgression(...
    'contourLevels', contourLevels, ...
    'intermediateIterationsToDisplay', [10 100], ...
    'displayedXrangeDegs', 0.25, ...
    'displayedYrangeDegs', 0.25 ...
    );

% Export to PDF
% NicePlot.exportFigToPDF('mosaicConstruction.pdf', hFig, 300);
    
% Show the final LMS mosaic
visualizedAperture = 'both'; % choose between 'both', 'lightCollectingArea', 'geometricArea'
theHexMosaic.visualizeGrid('visualizedConeAperture', visualizedAperture, 'generateNewFigure', true);
end


