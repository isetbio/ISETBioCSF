function run_stimulusAreaVaryConditions

    %% Inference engine and spatial summation sigma
    thresholdSignal = 'isomerizations';  % choose from {'isomerizations', 'photocurrents'}
    thresholdMethod = 'mlptGaussianRF';  % choose from {'mlpt', 'mlptGaussianRF'}
    
    
    %% Optics to employ. Choose from:
    % 'None': delta function PSF
    % 'AOoptics' : diffraction-limited with forced 6mm pupil
    % 'WvfHuman' : default human wavefront - based optics with mean (across subjects) Z-coeffs
    % 'WvfHumanMeanOTFmagMeanOTFphase' : human wavefront - based optics with mean (across subjects) OTF
    % 'Geisler' : Geisler optics
    % 'DavilaGeislerLsfAsPsf'
    % 'DavilaGeisler'
    
    employedOptics = 'AOoptics';

    spatialSummationData = containers.Map();
    spatialPoolingSigmaArcMinList = [1.5];
    for k = 1:numel(spatialPoolingSigmaArcMinList)
        spatialPoolingSigmaArcMin = spatialPoolingSigmaArcMinList(k);
        spatialSummationData(sprintf('summation_sigma_ArcMin_%2.2f',spatialPoolingSigmaArcMin)) = runSingleCondition(thresholdSignal, thresholdMethod, spatialPoolingSigmaArcMin, employedOptics);
    end
    
    save(sprintf('SummationDataExtendedRange_%s.mat',employedOptics), 'spatialSummationData', 'spatialPoolingSigmaArcMinList');
end

function plotData = runSingleCondition(thresholdSignal, thresholdMethod, spatialPoolingSigmaArcMin, employedOptics)
    
    fprintf('Running summation sigma: %2.3f arc min', spatialPoolingSigmaArcMin);
    [theDir,~] = fileparts(which(mfilename()));
    cd(theDir);
    
    %% Assemble all simulation params in a struct
    params = getParamsForStimulusAresVaryConditions(thresholdSignal, thresholdMethod, spatialPoolingSigmaArcMin, employedOptics);
    
    %% Simulation steps to perform
    % Re-compute the mosaic (true) or load it from the disk (false)
    params.computeMosaic = true; 
    
    % Re-compute the responses (true) or load them from the disk (false)
    params.computeResponses = true;
    
    % Compute photocurrents as well?
    params.computePhotocurrentResponseInstances = ~true;
    
    % Find the performance (true) or load performance data from the disk (false)
    params.findPerformance = true;
    
    % Fit the psychometric function? Set to true to obtain the threshols
    params.fitPsychometric = true;
    
    params.ramPercentageEmployed = 0.95;
    
    % Do not use the default plotting routine
    params.plotSpatialSummation = false;  % we will do our own plotting
    
    %% VISUALIZATION PARAMS
    visualizationScheme = {...
        %'mosaic' ...
        %'responses' ...
        'pooledSignal' ...
        'performance' ...
        'spatialScheme' ...    % graphic generated only during response computation
        'oiSequence' ...
        %'mosaic+emPath' ...    % graphic generated  only during response computation
        %'oiSequence' ...       % graphic generated  only during response computatio
    };
    params = visualizationParams(params, visualizationScheme);

    %% Go !
    [dataOut, extraData, rParams] = c_DavilaGeislerReplicateEyeMovements(params);
    
    %% Plot thresholds
    if strcmp(thresholdMethod,'mlpt')
        figurePDFname = sprintf('isomerizationsNoSpatialSummarion_Optics%s.pdf', employedOptics);
    else
        figurePDFname = sprintf('isomerizationsSpatialSummationSigma%2.2fArcMin_Optics%s.pdf',spatialPoolingSigmaArcMin, employedOptics);
    end
    
    plotData = generateThresholdPlot(dataOut, rParams, figurePDFname);
end

% ----- HELPER ROUTINES -----

function params = getParamsForStimulusAresVaryConditions(thresholdSignal, thresholdMethod, spatialPoolingSigmaArcMin, employedOptics)

    %% STIMULUS PARAMS
    % Varied stimulus area (spot diameter in arc min)
    params.spotDiametersMinutes = [0.3568 0.6181 1 2.5 5 10 20];
   
    % Stimulus background in degs
    params.backgroundSizeDegs = 55/60;
    
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
    params.imagePixels = 1000;
    
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
    
    %% OPTICS model and pupil size
    params.opticsModel = employedOptics;
    if (strcmpi(lower(params.opticsModel), 'aooptics'))
        params.pupilDiamMm = 6;
    else
        params.pupilDiamMm = 3;
    end
    
    %% MOSAIC PARAMS
    % Use a regularly-packed hegagonal mosaic (available packing options: 'rect', 'hex', 'hexReg')
    params.conePacking = 'hexReg';
    
    % Mosaic rotation in degrees
    params.mosaicRotationDegs = 0;
    
    % Cone spacing in microns
    params.coneSpacingMicrons = 3.0;
    
    % Cone aperture (inner-segment) diameter in microns
	params.innerSegmentSizeMicrons = 3.0;
    
    % Apply aperture low-pass
	params.apertureBlur = true;
    
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
	params.thresholdSignal = thresholdSignal;
    
    % Which inference engine to use to measure performance? Options are: 'mlpt', 'svm', 'svmSpaceTimeSeparable', 'svmGaussianRF', 'mlpt', 'mlpe'})); Threshold method to use
    params.thresholdMethod = thresholdMethod;

    % What performance (percent correct) do we require to declare that we reached the visibility threshold?
    params.thresholdCriterionFraction = 0.75;
    
    % If we are applying a PCA preprocessing step, how many spatiotemporal components to use?
    params.thresholdPCA = 60;
    
    %% Spatial pooling params
    params.spatialPoolingKernelParams = struct(...
        'subtractMeanOfNullResponseBeforeSummation', false, ...
        'type',  'GaussianRF', ...
        'shrinkageFactor', -spatialPoolingSigmaArcMin/60, ...  % If positive, sigma = spatialPoolingExtent * stimulus size, If negative, sigma =-spatialPoolingExtent
        'activationFunction', 'linear', ...
        'temporalPCAcoeffs', Inf, ...   ;  % Inf, results in no PCA, just the raw time series
        'adjustForConeDensity', false);
    
    
    % Fastrun  params
    %params.spotDiametersMinutes = [0.3568];
    %params.backgroundSizeDegs = 10/60;
    %params.nContrastsPerDirection = 2; 
    %params.nTrainingSamples = 2;
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

%% Visualize the performace curve
function plotData = generateThresholdPlot(dataOut, params, figurePDFname)
    % Plot data
    if strcmp(params.thresholdParams.method, 'mlpt')
        spatialPoolingSigmaMinArc = nan;
    else
        if (params.thresholdParams.spatialPoolingKernelParams.shrinkageFactor<0)
            spatialPoolingSigmaMinArc = abs(params.thresholdParams.spatialPoolingKernelParams.shrinkageFactor)*60;
        else
            spatialPoolingSigmaMinArc = params.thresholdParams.spatialPoolingKernelParams.shrinkageFactor * dataOut.spotDiametersMinutes;
        end
    end
    
    spotAreasMin2 = pi*((dataOut.spotDiametersMinutes/2).^2);
    summationAreaMin2 = pi * (2*spatialPoolingSigmaMinArc)^2;
    maxThresholdEnergies = params.maxSpotLuminanceCdM2 * params.temporalParams.stimulusDurationInSeconds*spotAreasMin2;

    lumIndex = 1;
    thresholdContrasts = [dataOut.mlptThresholds(lumIndex,:).thresholdContrasts];
    thresholdEnergies = thresholdContrasts.*maxThresholdEnergies;
    
    thresholdEnergyRange = [0.03 10];
    thresholdRange = [3*1e-5 1*1e-1];
    
    hFig = figure(100); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1000 450]);
    subplot(1,2,1);
    % Add Davila Geisler curve
    downShift = 0.8;
    A = LoadDigitizedDavilaGeislerFigure2;
    A(:,2) = (10^-downShift)*A(:,2);
    plot(A(:,1),A(:,2),'k.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.5 0.5 0.5],'LineWidth',1.5);
    hold on;
    plot(spotAreasMin2, thresholdEnergies, 'bo-', 'MarkerSize', 14, 'MarkerFaceColor', [0.7 0.9 1], 'LineWidth', 1.5);
    if (~isnan(spatialPoolingSigmaMinArc))
        plot(summationAreaMin2*[1 1], [thresholdEnergyRange(1) thresholdEnergyRange(2)], 'r-',  'LineWidth', 1.5);
    end
    
    axis 'square'; grid on;
    set(gca,'XScale','log','YScale','log', 'YLim', thresholdEnergyRange, 'FontSize', 16);
    xlabel('log10 spot area (square arc minutes)', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('log10 threshold energy', 'FontSize', 18, 'FontWeight', 'bold');
    legend({'davilaGeisler', 'isetbio', '$$A_{\mbox{sum}} = \pi(2\sigma_{\mbox{sum}})^2$'}, 'Interpreter', 'latex', 'Location', 'NorthWest');
    
    subplot(1,2,2);
    plot(spotAreasMin2, thresholdContrasts, 'bo-', 'MarkerSize', 14, 'MarkerFaceColor', [0.7 0.9 1], 'LineWidth', 1.5);
    hold on;
    plot(summationAreaMin2*[1 1],  [thresholdRange(1) thresholdRange(2)], 'r-',  'LineWidth', 1.5);
    axis 'square'; grid on;
    set(gca,'XScale','log','YScale','log',  'YLim', thresholdRange, 'FontSize', 16);
    xlabel('log10 spot area (square arc minutes)', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('log10 threshold contrast (arbitrary units)', 'FontSize', 18, 'FontWeight', 'bold');
    legend({'isetbio', '$A_{\mbox{sum}} = \pi(2\sigma_{\mbox{sum}})^2$'}, 'Interpreter', 'latex', 'Location', 'NorthWest');
    
    drawnow;
    NicePlot.exportFigToPDF(figurePDFname, hFig, 300);
    
    % Return data used in plot
    plotData.DavilaGeislerFigure.spotAras = squeeze(A(:,1));
    plotData.DavilaGeislerFigure.thresholdsEnergy = squeeze(A(:,2));
    plotData.spotArea = spotAreasMin2;
    plotData.thresholdsEnergy = thresholdEnergies;
    plotData.thresholdContrasts = thresholdContrasts;
    plotData.summationArea = summationAreaMin2;
end
