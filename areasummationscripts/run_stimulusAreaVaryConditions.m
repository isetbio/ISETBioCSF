function run_stimulusAreaVaryConditions

    %% Simulation params
    params = getParamsForStimulusAresVaryConditions();
    
    %% Simulation steps to perform
    % Re-compute the mosaic (true) or load it from the disk (false)
    params.computeMosaic = ~true; 
    
    % Re-compute the responses (true) or load them from the disk (false)
    params.computeResponses = ~true;
    
    % Compute photocurrents as well?
    params.computePhotocurrentResponseInstances = ~true;
    
    % Find the performance (true) or load performance data from the disk (false)
    params.findPerformance = true;
    
    % Fit the psychometric function? Set to true to obtain the threshols
    params.fitPsychometric = true;
    
    %% VISUALIZATION PARAMS
    visualizationScheme = {...
        %'mosaic' ...
        %'responses' ...
        'pooledSignal' ...
        %'performance' ...
        %'spatialScheme' ...    % graphic generated only during response computation
        %'mosaic+emPath' ...    % graphic generated  only during response computation
        %'oiSequence' ...       % graphic generated  only during response computation

    };
    params = visualizationParams(params, visualizationScheme);

    % Go
    [dataOut, extraData] = c_DavilaGeislerReplicateEyeMovements(params);
    dataOut
end


function params = getParamsForStimulusAresVaryConditions()

    %% STIMULUS PARAMS
    % Varied stimulus area (spot diameter in arc min)
    params.spotDiametersMinutes = [0.5 1 5 10 20];
    
    % Stimulus background in degs
    params.backgroundSizeDegs = 30/60;
    
    % Stimulus wavelength in nm
    params.wavelength = 550;
    
    % Stimulus luminance in cd/m2
    params.luminances = [10];
    
    % Stimulus duration in milliseconds
    params.stimDurationMs = 100;
    
    % Stimulus refresh rate in Ha
    params.stimRefreshRateHz = 50;
    
    % Add simulation time before stimulus onset (to stabilize the photocurrent)
	params.responseStabilizationMilliseconds = 10;
    
    % Add simulation time after stimulus offset (to let the photocurrent return to baserate)
	params.responseExtinctionMilliseconds = 20;

    % How many pixels to use to same the stimulus
    params.imagePixels = 400;
    
    % Lowest examined stimulus contrast
    params.lowContrast = 1e-6;
    
    % Highest examined stimulus contrast
    params.highContrast = 1e-2;
    
    % How many contrasts to use for the psychometric curve
    params.nContrastsPerDirection = 30;
    
	% How to space the constrasts, linearly or logarithmically?
    params.contrastScale = 'log';

    % Response instances to generate
    params.nTrainingSamples = 1024;
    
    
    %% OPTICS PARAMS
    % Pupil diameter in mm
    params.pupilDiamMm = 3;
    
    % Apply default human optics ?
    params.blur = true;
	params.opticsModel = 'WvfHuman';
    
    
    %% MOSAIC PARAMS
    % Use a regularly-packed hegagonal mosaic (available packing options: 'rect', 'hex', 'hexReg')
    params.conePacking = 'hexReg';
    
    % Mosaic rotation in degrees
    params.mosaicRotationDegs = 0;
    
    % Cone spacing in microns
    params.coneSpacingMicrons = 3.0;
    
    % Cone aperture (inner-segment) diameter in microns
	params.innerSegmentSizeMicrons = 3.0;
    
    % Apply aperture low-pass filtering ?
	params.apertureBlur = false;
    
    % Spatial density of L, M and S cones (sum: 1.0)
    params.LMSRatio = [0.67 0.33 0];

    % S-cone specific params
    % min distance between neighboring S-cones = f * local cone separation - to make the S-cone lattice semi-regular
    params.sConeMinDistanceFactor =  3.0;
    % radius of the s-cone free region (in microns)
	params.sConeFreeRadiusMicrons = 45;

    % Integration time (in seconds). Also the temporal support for the computed responses.
    params.integrationTime = 5.0/1000;
    
    % Dark noise, anyone?
    params.coneDarkNoiseRate = [0 0 0];

    % Freze the random noise generator so as to obtain repeatable results?
    params.freezeNoise = true;
    
    
    %% Eye movement params
    % What type of eye movement path to use (options are: 'frozen','frozen0', or 'random')
    params.emPathType = 'frozen0';
    
    % Center the eye movement path around (0,0)?
	params.centeredEMPaths = false;
    
    %% Performance params
    % What signal to use to measure performance? Options are: 'isomerizations', 'photocurrents'
	params.thresholdSignal = 'isomerizations';
    
    % Which inference engine to use to measure performance? Options are: 'mlpt', 'svm', 'svmSpaceTimeSeparable', 'svmGaussianRF', 'mlpt', 'mlpe'})); Threshold method to use
    params.thresholdMethod =  'mlpt';

    % What performance (percent correct) do we require to declare that we reached the visibility threshold?
    params.thresholdCriterionFraction = 0.75;
    
    % If we are applying a PCA preprocessing step, how many spatiotemporal components to use?
    params.thresholdPCA = 60;
    
    %% Spatial pooling params
    % default spatialPoolingKernel for use with 'svmGaussianRF' treshold method
    params.spatialPoolingKernelParams = struct(...
        'type',  'GaussianRF', ...
        'shrinkageFactor', 0.75', ...
        'activationFunction', 'linear', ...
        'temporalPCAcoeffs', Inf, ...   ;  % Inf, results in no PCA, just the raw time series
        'adjustForConeDensity', false);
end

function params = visualizationParams(params, visualizationScheme)

    availableVisualizations = {...
        'mosaic', ...
        'mosaic+emPath', ...
        'oiSequence', ...
        'responses', ...
        'spatialScheme', ...
        'pooledSignal', ...
        'performance' ...
    };
    

    for k = 1:numel(visualizationScheme)
        assert(ismember(visualizationScheme{k}, availableVisualizations), sprintf('Visualization ''%s'' is not available', visualizationScheme{k}));
    end
    
    % Visualize the mosaic ?
    params.visualizeMosaic = ismember('mosaic', visualizationScheme);
    % Visuzlize the mosaic with the first eye movement path superimposed?
    params.visualizeMosaicWithFirstEMpath =  ismember('mosaic+emPath', visualizationScheme);
    % Visualize the optical image sequence ?
    params.visualizeOIsequence = ismember('oiSequence', visualizationScheme);
    % Visualize the responses ?
    params.visualizeResponses = ismember('responses', visualizationScheme);
    % Visualize the spatial scheme (mosaic + stimulus)?
    params.visualizeSpatialScheme = ismember('spatialScheme', visualizationScheme);
    % Visualize the kernel transformed (spatial pooling) signals?
    params.visualizeKernelTransformedSignals = ismember('pooledSignal', visualizationScheme);
    % Visualize the performance ?
    params.visualizePerformance = ismember('performance', visualizationScheme);
end