function [validationData, extraData] = c_DavilaGeislerReplicateEyeMovements(varargin)
% c_DavilaGeislerReplicate(varargin)
%
% Compute thresholds to replicate spatial summation calculations of Davila and Geisler, more or less.
%
% This looks at thresholds as a function of spot size.  Our stimuli are
% monochromatic rather than monitor based, but to first order that should not make much difference
%
% Key/value pairs
%   'useScratchTopLevelDirName'- true/false (default false). 
%      When true, the top level output directory is [scratch]. 
%      When false, it is the name of this script.
%   'nTrainingSamples' - value (default 500).  Number of training samples to cycle through.
%   'spotDiametersMinutes' - vector (default [0.5 1 5 10 20 40]). Diameters of spots to be investigated.
%   'backgroundSizeDegs' - value (default 2.1). Size of square background
%   'wavelength' - value (default 550). Wavelength to use in calculations
%   'luminances' - vector (default [10]).  Background luminances in cd/m2 to be investigated.
%   'pupilDiamMm' - value (default 3).  Pupil diameter in mm.
%   'durationMs' - value (default 100).  Stimulus duration in msec.
%   'blur' - true/false (default true). Incorporate lens blur.
%   'opticsModel - string.  What optics model to use
%       'WvfHuman'           Isetbio standard wavefront based model of human optics (default).
%       'DavilaGeisler'      PSF based on DavilaGeisler line spread function.
%       'DavilaGeislerLsfAsPsf' Take D/G lsf and treat it directly as a psf
%       'Westheimer'         PSF based on Westheimer line spread function.
%       'Williams'           PSF based on Williams et al. MTF  
%   'innerSegmentSizeMicrons' - Diameter of the cone light-collecting area, in microns 
%       Default: 3.0, where 3 microns = 6 min arc for 300 mirons/degree in the human retina.
%   'apertureBlur' - Blur by cone aperture? true/false (default false).
%   'coneSpacingMicrons' - Cone spacing in microns (3).
%   'mosaicRotationDegs' - Rotation of hex or hexReg mosaic in degrees (default 0).
%   'coneDarkNoiseRate' - Vector of LMS dark (thermal) isomerization rate iso/sec (default [0,0,0]).
%   'LMSRatio' - Ratio of LMS cones in mosaic.  Should sum to 1 (default [0.67 0.33 0]).
%   'conePacking'   - how cones are packed spatially. 
%       'rect'              Rectangular mosaic.
%       'hex'               Hex mosaic with an eccentricity-varying cone spacing
%       'hexReg'            Hex mosaic witha regular cone spacing (default).
%   'imagePixels' - value (default 400).  Size of image pixel array
%   'nContrastsPerDirection' - value (default 30). Number of contrasts.
%   'lowContrast' - value (default 0.0001). Low contrast.
%   'highContrast' - value (default 0.1). High contrast.
%   'contrastScale' - 'log'/'linear' (default 'log'). Contrast scale.
%   'computeResponses' - true/false (default true).  Compute responses.
%   'findPerformance' - true/false (default true).  Find performance.
%   'thresholdMethod' - string (default 'mlpt').  How to find performance ('mlpt', 'mlpe', 'svm')
%   'thresholdPCA' - value (default 60).  Number of PCA components to keep for SVM method.
%   'fitPsychometric' - true/false (default true).  Fit psychometric functions.
%   'thresholdCriterionFraction' value (default 0.75). Criterion corrrect for threshold.
%   'generatePlots' - true/false (default true).  No plots are generated unless this is true.
%   'visualizeResponses' - true/false (default false).  Do the fancy response visualization when generating responses.
%   'visualizedResponseNormalization' - string (default 'submosaicBasedZscore'). How to normalize visualized responses
%        'submosaicBasedZscore'       [DHB Note:] Say what these do, please (default).
%        'LMSabsoluteResponseBased'
%        'LMabsoluteResponseBased'
%        'MabsoluteResponseBased'
%   'plotPsychometric' - true/false (default true).  Plot psychometric functions.
%   'plotSpatialSummation' - true/false (default true).  Plot results.
%   'freezeNoise' - true/false (default true). Freeze noise so calculations reproduce.

%% Parse input
p = inputParser;

%% -- PARAMS DEFINED IN c_DavilaGeislerReplicate (as of April 27, 2017) ---
p.addParameter('useScratchTopLevelDirName', false, @islogical);
p.addParameter('nTrainingSamples',500,@isnumeric);
p.addParameter('spotDiametersMinutes',[0.5 1 5 10 20 40],@isnumeric);
p.addParameter('backgroundSizeDegs',2.1,@isnumeric);
p.addParameter('wavelength',550,@isnumeric);
p.addParameter('luminances',[10],@isnumeric);
p.addParameter('pupilDiamMm',3,@isnumeric);
p.addParameter('durationMs',100,@isnumeric);
p.addParameter('blur',true,@islogical);
p.addParameter('opticsModel','WvfHuman',@ischar);
p.addParameter('innerSegmentSizeMicrons',3.0, @isnumeric);   % 3 microns = 0.6 min arc for 300 microns/deg in human retina
p.addParameter('apertureBlur', false, @islogical);
p.addParameter('coneSpacingMicrons', 3.0, @isnumeric);
p.addParameter('mosaicRotationDegs', 0, @isnumeric);
p.addParameter('coneDarkNoiseRate',[0 0 0], @isnumeric);
p.addParameter('LMSRatio',[0.67 0.33 0],@isnumeric);
p.addParameter('conePacking', 'hexReg',@ischar);                 
p.addParameter('imagePixels',400,@isnumeric);
p.addParameter('nContrastsPerDirection',30,@isnumeric);
p.addParameter('lowContrast',1e-6,@isnumeric);
p.addParameter('highContrast',1e-2,@isnumeric);
p.addParameter('contrastScale','log',@ischar);
p.addParameter('computeResponses',true,@islogical);
p.addParameter('findPerformance',true,@islogical);
%p.addParameter('thresholdMethod','mlpt',@ischar);
p.addParameter('thresholdPCA',60,@isnumeric);
p.addParameter('fitPsychometric',true,@islogical);
p.addParameter('thresholdCriterionFraction',0.75,@isnumeric);
p.addParameter('generatePlots',true,@islogical);
p.addParameter('visualizeResponses',false,@islogical);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.addParameter('plotPsychometric',true,@islogical);
p.addParameter('plotSpatialSummation',true,@islogical);
p.addParameter('freezeNoise',true,@islogical);

%% --- additional params for eye movements -----
% RESPONSE COMPUTATION OPTIONS
p.addParameter('emPathType','frozen0',@(x)ismember(x, {'none', 'frozen', 'frozen0', 'random'}));
p.addParameter('centeredEMPaths',false, @islogical); 
p.addParameter('responseStabilizationMilliseconds', 80, @isnumeric);
p.addParameter('responseExtinctionMilliseconds', 200, @isnumeric);
p.addParameter('computeMosaic',false,@islogical);

% MOSAIC OPTIONS
p.addParameter('integrationTime', 5.0/1000, @isnumeric);

% HEX MOSAIC OPTIONS
p.addParameter('sConeMinDistanceFactor', 3.0, @isnumeric); % min distance between neighboring S-cones = f * local cone separation - to make the S-cone lattice semi-regular
p.addParameter('sConeFreeRadiusMicrons', 45, @isnumeric);
p.addParameter('latticeAdjustmentPositionalToleranceF', 0.01, @isnumeric);
p.addParameter('latticeAdjustmentDelaunayToleranceF', 0.001, @isnumeric);
p.addParameter('marginF', [], @isnumeric);     

% DIAGNOSTIC OPTIONS
p.addParameter('displayTrialBlockPartitionDiagnostics', true, @islogical);
p.addParameter('displayResponseComputationProgress', false, @islogical);

% RESPONSE MAP VISUALIZATION OPTIONS
p.addParameter('visualizeMosaic',true, @islogical); 
p.addParameter('visualizeSpatialScheme', false, @islogical);
p.addParameter('visualizeKernelTransformedSignals',false, @islogical);
p.addParameter('visualizationFormat', 'montage', @ischar);
p.addParameter('visualizePerformance', false, @islogical);

% PERFORMANCE COMPUTATION OPTIONS
% default spatialPoolingKernel for use with 'svmGaussianRF' treshold method
spatialPoolingKernelParams = struct(...
    'type',  'GaussianRF', ...
    'shrinkageFactor', 0.75', ...
    'activationFunction', 'linear', ...
    'temporalPCAcoeffs', Inf, ...   ;  % Inf, results in no PCA, just the raw time series
    'adjustForConeDensity', false);
    
p.addParameter('spatialPoolingKernelParams', spatialPoolingKernelParams, @isstruct);
p.addParameter('useRBFSVMKernel', false, @islogical);
p.addParameter('thresholdMethod', 'mlpt', @(x)ismember(x, {'svm', 'svmSpaceTimeSeparable', 'svmGaussianRF', 'mlpt', 'mlpe'}));
p.addParameter('thresholdSignal', 'isomerizations', @(x)ismember(x, {'isomerizations', 'photocurrents'}));

p.parse(varargin{:});

%% Get the parameters we need
%
% Start with default
rParams = responseParamsGenerate('spatialType','spot','backgroundType','AO','modulationType','AO');

%% Set the  topLevelDir name
if (~p.Results.useScratchTopLevelDirName)
    rParams.topLevelDirParams.name = mfilename;
end

%% Loop over spatial frequency
for ll = 1:length(p.Results.luminances)
    for dd = 1:length(p.Results.spotDiametersMinutes)
        
        % Get stimulus parameters correct
        %
        % Spatial 
        rParams.spatialParams.spotSizeDegs = p.Results.spotDiametersMinutes(dd)/60;
        rParams.spatialParams.backgroundSizeDegs = p.Results.backgroundSizeDegs;
        rParams.spatialParams.fieldOfViewDegs = 1.1*p.Results.backgroundSizeDegs;
        rParams.spatialParams.row = p.Results.imagePixels;
        rParams.spatialParams.col = p.Results.imagePixels;
        rParams.spatialParams.viewingDistance = 7.16;  
                         
        % Set background wavelength
        rParams.backgroundParams.backgroundWavelengthsNm = [p.Results.wavelength];
        
        % Scale radiance to produce desired background levels
        deltaWl = 10;
        S = [p.Results.wavelength deltaWl 2];
        wls = SToWls(S);
        load T_xyz1931
        T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
        radiancePerUW = AOMonochromaticCornealPowerToRadiance(wls,rParams.backgroundParams.backgroundWavelengthsNm,1,rParams.oiParams.pupilDiamMm,rParams.spatialParams.backgroundSizeDegs^2);
        xyzPerUW = T_xyz*radiancePerUW*(wls(2)-wls(1));
        desiredLumCdM2 = p.Results.luminances(ll);
        rParams.backgroundParams.backgroundCornealPowerUW = desiredLumCdM2/xyzPerUW(2);
        
        % Scale spot parameters into what we think a reasonable range is.
        % In Davila & Geisler, they give thresholds in what they call
        % threshold energy, which is spotLuminance*durationSecs*areaMinutes2.
        %
        % Just to get things into a reasonable scale, we'll find a
        % luminance corresponding to a threshold energy of 10 for our
        % smallest spot.
        maxThresholdEnergy = 10;
        minSpotAreaMin2 = pi*(min(p.Results.spotDiametersMinutes)/2)^2;
        maxSpotLuminanceCdM2 = maxThresholdEnergy/((p.Results.durationMs/1000)*minSpotAreaMin2);
        rParams.colorModulationParams.startWl = p.Results.wavelength;
        rParams.colorModulationParams.endWl = p.Results.wavelength+deltaWl;
        rParams.colorModulationParams.deltaWl = deltaWl;
        rParams.colorModulationParams.spotWavelengthNm = p.Results.wavelength;
        rParams.colorModulationParams.spotCornealPowerUW = maxSpotLuminanceCdM2/xyzPerUW(2);
        rParams.colorModulationParams.contrast = 1;
        
        % Blur
        rParams.oiParams = modifyStructParams(rParams.oiParams, ...
        	'blur', p.Results.blur, ...
            'pupilDiamMm', p.Results.pupilDiamMm, ...
            'opticsModel', p.Results.opticsModel);
        
        % Their stimulus intervals were 100 msec each.
        %
        % Equate stimulusSamplingIntervalInSeconds to stimulusDurationInSeconds to generate 1 time point only.
        % Their main calculation was without eye movements
        stimulusDurationInSeconds = 100/1000;
        
        %%  --- ADDITIONS FOR EYE-MOVEMENTS / PHOTOCURRENTS ---
        % Assume a frame rate (in Hz). This is only used to sample the 100 ms stimulus pulse
        frameRate = 50;
        windowTauInSeconds = nan; % square-wave
        stimulusSamplingIntervalInSeconds = 1/frameRate;
        
        % Pre-stimulus time: allow photocurrent response to stabilize
        responseStabilizationSeconds = ceil(p.Results.responseStabilizationMilliseconds/1000/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
        % Post-stimulus time: allow photocurrent response to return to baseline
        responseExtinctionSeconds = ceil(p.Results.responseExtinctionMilliseconds/1000/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
        secondsToInclude = responseStabilizationSeconds+stimulusDurationInSeconds+responseExtinctionSeconds;
        rParams.temporalParams = modifyStructParams(rParams.temporalParams, ...
            'frameRate', frameRate, ...
            'windowTauInSeconds', windowTauInSeconds, ...
            'stimulusSamplingIntervalInSeconds', stimulusSamplingIntervalInSeconds, ...
            'stimulusDurationInSeconds', stimulusDurationInSeconds, ...
            'secondsToInclude', secondsToInclude, ...
            'secondsForResponseStabilization', responseStabilizationSeconds, ...
            'secondsForResponseExtinction', responseExtinctionSeconds, ...
            'secondsToIncludeOffset', 0/1000, ...
            'emPathType', p.Results.emPathType ...
        );
        %%  --- END OF ADDITIONS FOR EYE-MOVEMENTS / PHOTOCURRENTS ---
        
        
        % Set up mosaic parameters. NOTE THAT IN THIS VERSION (WITH EYE MOVEMENTS)
        % WE DO NOT INTEGRATE FOR THE ENTIRE STIMULUS 
        
        if (isfield(rParams.mosaicParams, 'realisticSconeSubmosaic'))
            error('The realisticSconeSubmosaic has been removed. Pass directly sConeMinDistanceFactor and sConeFreeRadiusMicrons\n');
        end
        
        rParams.mosaicParams = modifyStructParams(rParams.mosaicParams, ...
            'fieldOfViewDegs', 1.1*p.Results.backgroundSizeDegs, ...  
            'innerSegmentSizeMicrons',p.Results.innerSegmentSizeMicrons, ...
            'apertureBlur',p.Results.apertureBlur, ...
            'mosaicRotationDegs',p.Results.mosaicRotationDegs,...
            'coneDarkNoiseRate',p.Results.coneDarkNoiseRate,...
            'LMSRatio',p.Results.LMSRatio,...
            'coneSpacingMicrons', p.Results.coneSpacingMicrons, ...
            'conePacking', p.Results.conePacking, ...
        	'integrationTimeInSeconds', p.Results.integrationTime, ...   % <<<<<<<<< DIFFERENT FOR EYE MOVEMENTS
        	'isomerizationNoise', 'random',...              % select from {'random', 'frozen', 'none'}
        	'osNoise', 'random', ...                        % select from {'random', 'frozen', 'none'}
        	'osModel', 'Linear');
        
        if strcmp(p.Results.conePacking, 'hex')
            rParams.mosaicParams = modifyStructParams(rParams.mosaicParams, ...
                'sConeMinDistanceFactor', p.Results.sConeMinDistanceFactor, ...
                'sConeFreeRadiusMicrons', p.Results.sConeFreeRadiusMicrons, ...
                'latticeAdjustmentPositionalToleranceF', p.Results.latticeAdjustmentPositionalToleranceF, ...
                'latticeAdjustmentDelaunayToleranceF', p.Results.latticeAdjustmentDelaunayToleranceF, ...
                'marginF', p.Results.marginF);
        end
        
        % Make sure mosaic noise parameters match the frozen noise flag
        % passed in.  Doing this here means that things work on down the
        % chain, while if we don't then there is a problem because routines
        % create the wrong directory.
        if (p.Results.freezeNoise)
            if (strcmp(rParams.mosaicParams.isomerizationNoise, 'random'))
                rParams.mosaicParams.isomerizationNoise = 'frozen';
            end
            if (strcmp(rParams.mosaicParams.osNoise, 'random'))
                rParams.mosaicParams.osNoise = 'frozen';
            end
        end

        % Parameters that define the LM instances we'll generate here
        %
        % Use default LMPlane.
        testDirectionParams = instanceParamsGenerate('instanceType','contrasts');
        testDirectionParams = modifyStructParams(testDirectionParams, ...
            'trialsNum', p.Results.nTrainingSamples, ...
            'nContrastsPerDirection', p.Results.nContrastsPerDirection, ...
            'lowContrast', p.Results.lowContrast, ...
        	'highContrast', p.Results.highContrast, ...
        	'contrastScale', p.Results.contrastScale ...                       % choose between 'linear' and 'log'
            );
        
        % Parameters related to how we find thresholds from responses
        thresholdParams = thresholdParamsGenerate;

        %%  --- MODIFICATIONS FOR SVM SPATIAL SUMMATOR PREPROCESSING ---
        thresholdParams = modifyStructParams(thresholdParams, ...
        	'criterionFraction', p.Results.thresholdCriterionFraction, ...
        	'method', p.Results.thresholdMethod, ...
            'STANDARDIZE', false, ...                    % Standardize data for PCA
            'standardizeSVMpredictors', false, ...       % Standardize data for SVM
            'useRBFKernel', p.Results.useRBFSVMKernel, ...
        	'PCAComponents', p.Results.thresholdPCA, ...
            'signalSource', p.Results.thresholdSignal ...
            );
        if (strcmp(thresholdParams.method, 'svmGaussianRF'))
            thresholdParams = modifyStructParams(thresholdParams, ...
                'spatialPoolingKernelParams', p.Results.spatialPoolingKernelParams);  
        end
        
        
        %% Compute response instances
        if (p.Results.computeResponses)
           t_coneCurrentEyeMovementsResponseInstances(...
               'rParams',rParams,...
               'testDirectionParams',testDirectionParams,...
               'compute',true,...
               'computeMosaic', (ll == 1) && (dd == 1) && (p.Results.computeMosaic), ...
               'visualizedResponseNormalization', p.Results.visualizedResponseNormalization, ...
               'visualizeResponses',p.Results.visualizeResponses, ...
               'visualizeMosaic', p.Results.visualizeMosaic, ...
               'generatePlots',p.Results.generatePlots,...
               'freezeNoise',p.Results.freezeNoise, ...
               'centeredEMPaths',p.Results.centeredEMPaths, ...                 % new
               'visualizeSpatialScheme', p.Results.visualizeSpatialScheme); ...   % new
        end
        
        if ((p.Results.visualizeResponses) || (p.Results.visualizeSpatialScheme)) ||  (p.Results.visualizeMosaic) 
            t_coneCurrentEyeMovementsResponseInstances(...
               'rParams',rParams,...
               'testDirectionParams',testDirectionParams,...
               'compute', false,...
               'computeMosaic', false, ...
               'visualizedResponseNormalization', p.Results.visualizedResponseNormalization, ...
               'visualizeResponses',p.Results.visualizeResponses, ...
               'visualizeMosaic', p.Results.visualizeMosaic, ...
               'generatePlots', true,...
               'freezeNoise',p.Results.freezeNoise, ...
               'centeredEMPaths',p.Results.centeredEMPaths, ...                                 
               'visualizeSpatialScheme', p.Results.visualizeSpatialScheme, ...
               'visualizeOuterSegmentFilters', false);
        end
    
        %% Find performance, template max likeli
        if (p.Results.findPerformance) || (p.Results.visualizePerformance)
            t_colorDetectFindPerformance(...
                'rParams',rParams, ...
                'testDirectionParams',testDirectionParams, ...
                'thresholdParams',thresholdParams, ...
                'compute',p.Results.findPerformance, ...
                'visualizeSpatialScheme', p.Results.visualizeSpatialScheme, ...
                'visualizeKernelTransformedSignals', p.Results.visualizeKernelTransformedSignals, ...
                'plotSvmBoundary',false, ...
                'plotPsychometric',false, ...
                'freezeNoise',p.Results.freezeNoise);
        end
        
        %% Fit psychometric functions
        if (p.Results.fitPsychometric)
            davilaGeislerReplicate.spotDiametersMinutes(ll,dd) = p.Results.spotDiametersMinutes(dd);
            davilaGeislerReplicate.mlptThresholds(ll,dd) = ...
                t_fitPsychometricFunctions(...
                    'rParams',rParams,'instanceParams',...
                    testDirectionParams,'thresholdParams',...
                    thresholdParams, ...
                    'generatePlots',p.Results.generatePlots && p.Results.plotPsychometric);
            %close all;
        end
    end
end
%% Write out the data

%
% This read and write does not distinguish the number of backgrounds
% studied, and so can screw up if the number of backgrounds changes 
% between a run that generated the data and one that read it back.  If
% everything else is in place, it is quick to regenerate the data by just
% doing the fit psychometric step, which is pretty quick.
if (p.Results.fitPsychometric) && ((p.Results.findPerformance) || (p.Results.visualizePerformance))
    fprintf('Writing performance data ... ');
    paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams, thresholdParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    rwObject.write('davilaGeislerReplicateEyeMovements',davilaGeislerReplicate,paramsList,writeProgram);
    fprintf('done\n');

    %% Read data back in
    fprintf('Reading performance data ...');
    % nameParams = rParams.spatialParams;
    % nameParams.spotSizeDegs = 0;
    % nameParams.backgroundSizeDegs = 0;
    % nameParams.fieldOfViewDegs = 0;
    % paramsList = {rParams.topLevelDirParams, nameParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
    paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams, thresholdParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    davilaGeislerReplicate = rwObject.read('davilaGeislerReplicateEyeMovements',paramsList,writeProgram);
    fprintf('done\n');



    %% Output validation data
    if (nargout > 0)
        validationData.spotDiametersMinutes = davilaGeislerReplicate.spotDiametersMinutes;
        validationData.mlptThresholds = davilaGeislerReplicate.mlptThresholds;
        validationData.luminances = p.Results.luminances;
        extraData.paramsList = paramsList;
        extraData.p.Results = p.Results;
    end


    %% Make a plot of estimated threshold versus training set size
    %
    % The way the plot is coded counts on the test contrasts never changing
    % across the conditions, which we could explicitly check for here.
    if (p.Results.generatePlots && p.Results.plotSpatialSummation)
        % Some numbers we need
        spotAreasMin2 = pi*((p.Results.spotDiametersMinutes/2).^2);
        maxThresholdEnergies = maxSpotLuminanceCdM2*rParams.temporalParams.stimulusDurationInSeconds*spotAreasMin2;

        % Davila-Geisler style plot
        hFig = figure; clf; hold on
        fontBump = 4;
        markerBump = -4;
        set(gcf,'Position',[100 100 450 650]);
        set(gca,'FontSize', 14+fontBump);
        theColors = ['r' 'g' 'b'];
        legendStr = cell(length(p.Results.luminances),1);
        for ll = 1:length(p.Results.luminances)
            theColorIndex = rem(ll,length(theColors)) + 1;
            plot(spotAreasMin2,[davilaGeislerReplicate.mlptThresholds(ll,:).thresholdContrasts].*maxThresholdEnergies, ...
                [theColors(theColorIndex) 'o-'],'MarkerSize',12+markerBump,'MarkerFaceColor',theColors(theColorIndex),'LineWidth',3);  
            legendStr{ll} = sprintf('%0.1f cd/m2',p.Results.luminances(ll));
        end

        % Add Davila Geisler curve
        downShift = 0.8;
        A = LoadDigitizedDavilaGeislerFigure2;
        A(:,2) = (10^-downShift)*A(:,2);
        plot(A(:,1),A(:,2),'k:','LineWidth',1);

        set(gca,'XScale','log');
        set(gca,'YScale','log');
        xlabel('Log10 Spot Area (square arc minutes)', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
        ylabel('Log10 Threshold Energy', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
        xlim([1e-2 1e4]); ylim([1e-3 1e3]);
        legend(legendStr,'Location','NorthWest','FontSize',rParams.plotParams.labelFontSize+fontBump);
        box off; grid on
        if (p.Results.blur)
            title(sprintf('Computational Observer CSF - w/ blur',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
            rwObject.write('davilaGeislerReplicateWithBlur',hFig,paramsList,writeProgram,'Type','figure');
        else
            title(sprintf('Computational Observer CSF - no blur',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
            rwObject.write('davilaGeislerReplicateNoBlur',hFig,paramsList,writeProgram,'Type','figure');
        end

        % Plot it the way we plot our data
        hFig1 = figure;
        set(gcf,'Position',[100 100 450 650]);
        set(gca,'FontSize', rParams.plotParams.axisFontSize+fontBump);
        legendStr = cell(length(p.Results.luminances),1);
        for ll = 1:length(p.Results.luminances)
            theColorIndex = rem(ll,length(theColors)) + 1;
            plot(spotAreasMin2,[davilaGeislerReplicate.mlptThresholds(ll,:).thresholdContrasts], ...
                [theColors(theColorIndex) 'o-'],'MarkerSize',rParams.plotParams.markerSize+markerBump,'MarkerFaceColor',theColors(theColorIndex),'LineWidth',rParams.plotParams.lineWidth);  
            legendStr{ll} = sprintf('%0.1f cd/m2',p.Results.luminances(ll));
        end

        set(gca,'XScale','log');
        set(gca,'YScale','log');
         xlabel('Log10 Spot Area (square arc minutes)', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
        ylabel('Log10 Threshold Contrast (arb. units)', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
        xlim([1e-2 1e4]);
        ylim([1e-3 1e3]);
        axis('square');
        legend(legendStr,'Location','NorthWest','FontSize',rParams.plotParams.labelFontSize+fontBump);
        box off; grid on
        if (p.Results.blur)
            title(sprintf('Computational Observer CSF - w/ blur',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
            rwObject.write('davilaGeislerReplicateWithBlur1',hFig1,paramsList,writeProgram,'Type','figure');
        else
            title(sprintf('Computational Observer CSF - no blur',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
            rwObject.write('davilaGeislerReplicateNoBlur1',hFig1,paramsList,writeProgram,'Type','figure');
        end

        % Save data for plots in convenient form
        plotData.spotAreasMin2 = spotAreasMin2;
        plotData.thresholdEnergy = [davilaGeislerReplicate.mlptThresholds(:,:).thresholdContrasts].*maxThresholdEnergies;
        plotData.thresholdContrasts = [davilaGeislerReplicate.mlptThresholds(:,:).thresholdContrasts];
        plotData.digitizedAreasMin2 = A(:,1)';
        plotData.digitizedEnergy = A(:,2)';
        if (p.Results.blur)
            rwObject.write('davilaGeislerReplicatePlotDataWithBlur',plotData,paramsList,writeProgram);

        else
            rwObject.write('davilaGeislerReplicatePlotDataNoBlur',plotData,paramsList,writeProgram);
        end
    end
end

end

