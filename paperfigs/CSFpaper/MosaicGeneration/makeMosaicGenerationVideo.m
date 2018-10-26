function makeMosaicGenerationVideo
% cd to wherver this script resides
[localDir,~] = fileparts(which(mfilename()));
cd(localDir)

% Set random seed to obtain replicable results
rng(1235);

params.fovDegs = [0.25 0.25]; % FOV in degrees ([width height], default: 0.25x0.25
makeNew = true;

if (makeNew)
    % Set coneMosaicHex - specific params
    params.resamplingFactor = 15;                            % 9 is good;how fine to sample the hex mosaic positions with an underlying rect mosaic                         
    params.eccBasedConeDensity = true;                      % if true, generate a mosaic where cone spacing varies with eccentricity (Curcio model)
    params.customLambda = [];                               % cone spacing in microns (only used with regular hex mosaics)
    params.rotationDegs = 0;                                % rotation of the mosaic, in degrees 0 = , 30 = (makes sense with regular hex mosaics)
    params.customInnerSegmentDiameter = [];                 % inner segment diameter, in microns (empty for isetbio default) 
    params.spatialDensity = [0 0.6 0.3 0.1];                % K/L/M/S cone densities
    params.sConeMinDistanceFactor = 3.0;                    % min distance between neighboring S-cones = f * local cone separation - to make the S-cone lattice semi-regular
    params.sConeFreeRadiusMicrons = 45;                     % radius of S-cone free retina, in microns
    params.latticeAdjustmentPositionalToleranceF = [];      % determines cone delta movement tolerance for terminating iterative adjustment - by default this is 0.01 (here setting it lower for faster, but less acurate mosaic generation)
    params.latticeAdjustmentDelaunayToleranceF = [];        % determines position tolerance for triggering another Delaunay triangularization - by default this is 0.001 (here setting it lower for faster, but less acurate mosaic generation)
    params.maxGridAdjustmentIterations = 3000;
    params.marginF = [];
    saveLatticeAdjustmentProgression = true;                % set to true, only if interested to see how the mosaic lattice is iteratively adjusted when eccBasedConeDensity is true               

    % Generate the mosaic
    theHexMosaic = coneMosaicHex(params.resamplingFactor, ...
        'fovDegs', params.fovDegs, ...
        'eccBasedConeDensity', params.eccBasedConeDensity, ...
        'rotationDegs', params.rotationDegs, ...
        'spatialDensity', params.spatialDensity, ...
        'sConeMinDistanceFactor', params.sConeMinDistanceFactor, ...    
        'sConeFreeRadiusMicrons', params.sConeFreeRadiusMicrons, ...    
        'customLambda', params.customLambda, ...
        'customInnerSegmentDiameter', params.customInnerSegmentDiameter, ...
        'latticeAdjustmentPositionalToleranctbUeF', params.latticeAdjustmentPositionalToleranceF, ...              
        'latticeAdjustmentDelaunayToleranceF', params.latticeAdjustmentDelaunayToleranceF, ...     
        'marginF', params.marginF, ...
        'maxGridIterations', 5000, ...
        'saveLatticeAdjustmentProgression', saveLatticeAdjustmentProgression ...  
    );
    % Save the mosaic
    save(sprintf('theHexMosaic%2.2fdegs.mat',max(params.fovDegs)), 'theHexMosaic', '-v7.3');
else
    fprintf('\nLoading mosaic ... ');
    % Load the mosaic
    load(sprintf('theHexMosaic%2.2fdegs.mat',max(params.fovDegs)));
    fprintf('Done\n');
end

hFig = figure(1); clf;
set(hFig, 'Position', [10 10 1100 1100], 'Color', [1 1 1]);
subplot('Position', [0.03 0.01 0.94 0.94]);

lastIter = 890;
iterations = [0:50 50:10:100 100:50:500 550:100:lastIter lastIter];
radius = 1e-6*100;
videoOBJ = VideoWriter(sprintf('MosaicConstruction.mp4'), 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 4;
    videoOBJ.Quality = 100;
    videoOBJ.open();
    
for i = 1:numel(iterations)
    iteration = iterations(i);
    cla(gca);
    theHexMosaic.plotMosaicAtIteration(gca, iteration, radius);
    videoOBJ.writeVideo(getframe(hFig));
    if (iteration == 1)
        for k = 1:5
            videoOBJ.writeVideo(getframe(hFig));
        end
    end
end
videoOBJ.close();

end




