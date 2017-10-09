function varargout =  c_DavilaGeislerReplicateEyeMovements(varargin)
% c_DavilaGeislerReplicate(varargin)
%
% Compute thresholds to replicate spatial summation calculations of Davila and Geisler, more or less.
%
% This looks at thresholds as a function of spot size.  Our stimuli are
% monochromatic rather than monitor based, but to first order that should not make much difference

%% Set the default output
varargout = {};
varargout{1} = [];  % thresholds for different stimulus sizes and background luminances
varargout{2} = [];  % extraData

%% Parse the input params
userParams = parseUserParams(varargin{:});

% Generate runParams
[rParams, testDirectionParams, thresholdParams] = runParams(userParams);

%% Loop over spatial frequency
for backgroundLumIndex = 1:length(userParams.luminances)
    for stimDiamIndex = 1:length(userParams.spotDiametersMinutes)
        
        % Update rParams for current stimDiamIndex and backgroundLumIndex
        rParams = updateRunParams(rParams, userParams, backgroundLumIndex, stimDiamIndex);
 
        % Compute response instances
        if ((userParams.computeResponses) || ...
            (userParams.visualizeResponses) || (userParams.visualizeSpatialScheme) ||  (userParams.visualizeMosaic) ...     
           )
        % Call the core response computation routine to either generate or visualize the responses
        t_coneCurrentEyeMovementsResponseInstances(...
               'rParams',rParams,...
               'testDirectionParams',testDirectionParams,...
               'compute',userParams.computeResponses,...
               'computePhotocurrentResponseInstances', userParams.computePhotocurrentResponseInstances, ...
               'computeMosaic', (backgroundLumIndex == 1) && (stimDiamIndex == 1) && (userParams.computeMosaic), ...
               'freezeNoise', userParams.freezeNoise, ...
               'centeredEMPaths', userParams.centeredEMPaths, ... 
               'generatePlots', userParams.generatePlots,...
               'visualizeOIsequence', userParams.visualizeOIsequence, ....
               'visualizeMosaicWithFirstEMpath', userParams.visualizeMosaicWithFirstEMpath, ...
               'visualizedResponseNormalization', userParams.visualizedResponseNormalization, ...
               'visualizeResponses',userParams.visualizeResponses, ...
               'visualizeMosaic', userParams.visualizeMosaic, ...         
               'visualizeSpatialScheme', userParams.visualizeSpatialScheme,...
               'visualizeOuterSegmentFilters', false...
           ); 
        end
    
        % Find performance
        if (userParams.findPerformance) || (userParams.visualizePerformance)
            % Compute psychometric function
            t_colorDetectFindPerformance(...
                'rParams',rParams, ...
                'testDirectionParams',testDirectionParams, ...
                'thresholdParams',thresholdParams, ...
                'compute',userParams.findPerformance, ...
                'visualizeSpatialScheme', userParams.visualizeSpatialScheme, ...
                'visualizeKernelTransformedSignals', userParams.visualizeKernelTransformedSignals, ...
                'plotSvmBoundary',false, ...
                'plotPsychometric', (userParams.visualizePerformance) & (~userParams.fitPsychometric), ...
                'freezeNoise',userParams.freezeNoise);
        
            % Fit psychometric function
            if (userParams.fitPsychometric)
                davilaGeislerReplicate.spotDiametersMinutes(backgroundLumIndex,stimDiamIndex) = userParams.spotDiametersMinutes(stimDiamIndex);
                davilaGeislerReplicate.mlptThresholds(backgroundLumIndex,stimDiamIndex) = ...
                    t_fitPsychometricFunctions(...
                        'rParams',rParams,'instanceParams',...
                        testDirectionParams,'thresholdParams',...
                        thresholdParams, ...
                        'generatePlots',userParams.generatePlots && userParams.plotPsychometric);
                %close all;
            end
        end
    end % stimDiamIndex
end % backgroundLumIndex

[dataOut, extraData] = exportPerformanceData(rParams, userParams, thresholdParams, davilaGeislerReplicate, mfilename);
if (nargout > 0)
    varargout{1} = dataOut;
    varargout{2} = extraData;
end
end


% -------- HELPER FUNCTIONS ---------

function [dataOut, extraData] = exportPerformanceData(rParams, userParams, thresholdParams, davilaGeislerReplicate, writeProgram)
dataOut = [];
extraData = [];
    
% This read and write does not distinguish the number of backgrounds
% studied, and so can screw up if the number of backgrounds changes 
% between a run that generated the data and one that read it back.  If
% everything else is in place, it is quick to regenerate the data by just
% doing the fit psychometric step, which is pretty quick.
if (userParams.fitPsychometric) && ((userParams.findPerformance) || (userParams.visualizePerformance))
    fprintf('Writing performance data ... ');
    paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams, thresholdParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    rwObject.write('davilaGeislerReplicateEyeMovements',davilaGeislerReplicate,paramsList,writeProgram);
    fprintf('done\n');

    % Read data back in
    fprintf('Reading performance data ...');
    paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams, thresholdParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    davilaGeislerReplicate = rwObject.read('davilaGeislerReplicateEyeMovements',paramsList,writeProgram);
    fprintf('done\n');

    % Returned data
    dataOut.spotDiametersMinutes = davilaGeislerReplicate.spotDiametersMinutes;
    dataOut.mlptThresholds = davilaGeislerReplicate.mlptThresholds;
    dataOut.luminances = userParams.luminances;
    extraData.paramsList = paramsList;
    extraData.userParams = userParams;

    % Make a plot of estimated threshold versus training set size
    %
    % The way the plot is coded counts on the test contrasts never changing
    % across the conditions, which we could explicitly check for here.
    if (userParams.generatePlots && userParams.plotSpatialSummation)
        % Some numbers we need
        spotAreasMin2 = pi*((userParams.spotDiametersMinutes/2).^2);
        maxThresholdEnergies = rParams.maxSpotLuminanceCdM2*rParams.temporalParams.stimulusDurationInSeconds*spotAreasMin2;

        % Davila-Geisler style plot
        hFig = figure; clf; hold on
        fontBump = 4;
        markerBump = -4;
        set(gcf,'Position',[100 100 450 650]);
        set(gca,'FontSize', 14+fontBump);
        theColors = ['r' 'g' 'b'];
        legendStr = cell(length(userParams.luminances),1);
        for ll = 1:length(userParams.luminances)
            theColorIndex = rem(ll,length(theColors)) + 1;
            plot(spotAreasMin2,[davilaGeislerReplicate.mlptThresholds(ll,:).thresholdContrasts].*maxThresholdEnergies, ...
                [theColors(theColorIndex) 'o-'],'MarkerSize',12+markerBump,'MarkerFaceColor',theColors(theColorIndex),'LineWidth',3);  
            legendStr{ll} = sprintf('%0.1f cd/m2',userParams.luminances(ll));
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
        if (userParams.blur)
            title(sprintf('Computational Observer CSF - w/ blur'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
            rwObject.write('davilaGeislerReplicateWithBlur',hFig,paramsList,writeProgram,'Type','figure');
        else
            title(sprintf('Computational Observer CSF - no blur'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
            rwObject.write('davilaGeislerReplicateNoBlur',hFig,paramsList,writeProgram,'Type','figure');
        end

        % Plot it the way we plot our data
        hFig1 = figure;
        set(gcf,'Position',[100 100 450 650]);
        set(gca,'FontSize', rParams.plotParams.axisFontSize+fontBump);
        legendStr = cell(length(userParams.luminances),1);
        for ll = 1:length(userParams.luminances)
            theColorIndex = rem(ll,length(theColors)) + 1;
            plot(spotAreasMin2,[davilaGeislerReplicate.mlptThresholds(ll,:).thresholdContrasts], ...
                [theColors(theColorIndex) 'o-'],'MarkerSize',rParams.plotParams.markerSize+markerBump,'MarkerFaceColor',theColors(theColorIndex),'LineWidth',rParams.plotParams.lineWidth);  
            legendStr{ll} = sprintf('%0.1f cd/m2',userParams.luminances(ll));
        end

        set(gca,'XScale','log');
        set(gca,'YScale','log');
         xlabel('Log10 Spot Area (square arc minutes)', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
        ylabel('Log10 Threshold Contrast (arb. units)', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
        xlim([1e-2 1e4]);
        ylim([1e-5 1e1]);
        axis('square');
        legend(legendStr,'Location','NorthWest','FontSize',rParams.plotParams.labelFontSize+fontBump);
        box off; grid on
        if (userParams.blur)
            title(sprintf('Computational Observer CSF - w/ blur'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
            rwObject.write('davilaGeislerReplicateWithBlur1',hFig1,paramsList,writeProgram,'Type','figure');
        else
            title(sprintf('Computational Observer CSF - no blur'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
            rwObject.write('davilaGeislerReplicateNoBlur1',hFig1,paramsList,writeProgram,'Type','figure');
        end

        % Save data for plots in convenient form
        plotData.spotAreasMin2 = spotAreasMin2;
        plotData.thresholdEnergy = [davilaGeislerReplicate.mlptThresholds(:,:).thresholdContrasts].*maxThresholdEnergies;
        plotData.thresholdContrasts = [davilaGeislerReplicate.mlptThresholds(:,:).thresholdContrasts];
        plotData.digitizedAreasMin2 = A(:,1)';
        plotData.digitizedEnergy = A(:,2)';
        if (userParams.blur)
            rwObject.write('davilaGeislerReplicatePlotDataWithBlur',plotData,paramsList,writeProgram);

        else
            rwObject.write('davilaGeislerReplicatePlotDataNoBlur',plotData,paramsList,writeProgram);
        end
    end
end
end


function [rParams, testDirectionParams, thresholdParams] = runParams(userParams)
% Start with default
rParams = responseParamsGenerate('spatialType','spot','backgroundType','AO','modulationType','AO');

%% Set the  topLevelDir name
if (~userParams.useScratchTopLevelDirName)
    rParams.topLevelDirParams.name = mfilename;
end

%% Parameters that define the LM instances we'll generate here
testDirectionParams = instanceParamsGenerate('instanceType','contrasts');
testDirectionParams = modifyStructParams(testDirectionParams, ...
    'trialsNum', userParams.nTrainingSamples, ...
    'nContrastsPerDirection', userParams.nContrastsPerDirection, ...
    'lowContrast', userParams.lowContrast, ...
    'highContrast', userParams.highContrast, ...
    'contrastScale', userParams.contrastScale ...  
    );
    
%% Parameters related to how we find thresholds from responses
thresholdParams = thresholdParamsGenerate;

%%  --- MODIFICATIONS FOR SVM SPATIAL SUMMATOR PREPROCESSING ---
thresholdParams = modifyStructParams(thresholdParams, ...
    'criterionFraction', userParams.thresholdCriterionFraction, ...
    'method', userParams.thresholdMethod, ...
    'spatialPoolingKernelParams', userParams.spatialPoolingKernelParams, ...
    'STANDARDIZE', false, ...                    % Standardize data for PCA
    'standardizeSVMpredictors', false, ...       % Standardize data for SVM
    'useRBFKernel', userParams.useRBFSVMKernel, ...
    'PCAComponents', userParams.thresholdPCA, ...
    'signalSource', userParams.thresholdSignal ...
    )
userParams.spatialPoolingKernelParams
if (strcmp(thresholdParams.method, 'svmGaussianRF')) 
    thresholdParams = modifyStructParams(thresholdParams, ...
        'spatialPoolingKernelParams', userParams.spatialPoolingKernelParams);  
end       
end



function userParams = parseUserParams(paramsList)

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
%   'computePhotocurrentResponseInstances' - true/false (default false). Compute photocurrent responses
%   'findPerformance' - true/false (default true).  Find performance.
%   'thresholdMethod' - string (default 'mlpt').  How to find performance ('mlpt', 'mlpe', 'svm')
%   'thresholdPCA' - value (default 60).  Number of PCA components to keep for SVM method.
%   'fitPsychometric' - true/false (default true).  Fit psychometric functions.
%   'thresholdCriterionFraction' value (default 0.75). Criterion corrrect for threshold.
%   'generatePlots' - true/false (default true).  No plots are generated unless this is true.
%   'visualizeOIsequence' - true/false (default false). Visualize the sequence of optical images
%   'visualizeMosaicWithFirstEMpath' - true/false (default false). Visualize the mosaic with the first EMpath superimposed
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
p.addParameter('nTrainingSamples',1024,@isnumeric);
p.addParameter('spotDiametersMinutes',[0.5 1 5 10 20],@isnumeric);
p.addParameter('backgroundSizeDegs',30/60,@isnumeric);
p.addParameter('wavelength',550,@isnumeric);
p.addParameter('luminances',[10],@isnumeric);
p.addParameter('pupilDiamMm',3,@isnumeric);
p.addParameter('stimDurationMs',100,@isnumeric);
p.addParameter('stimRefreshRateHz', 50, @isnumeric);

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
p.addParameter('computePhotocurrentResponseInstances', false, @islogical);
p.addParameter('findPerformance',true,@islogical);
p.addParameter('thresholdPCA',60,@isnumeric);
p.addParameter('fitPsychometric',true,@islogical);
p.addParameter('thresholdCriterionFraction',0.75,@isnumeric);
p.addParameter('generatePlots',true,@islogical);
p.addParameter('visualizeResponses',true,@islogical);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.addParameter('plotPsychometric',true,@islogical);
p.addParameter('plotSpatialSummation',true,@islogical);
p.addParameter('freezeNoise',true,@islogical);

%% --- additional params for eye movements -----
% RESPONSE COMPUTATION OPTIONS
p.addParameter('emPathType','frozen0',@(x)ismember(x, {'none', 'frozen', 'frozen0', 'random'}));
p.addParameter('centeredEMPaths',false, @islogical); 
p.addParameter('responseStabilizationMilliseconds', 10, @isnumeric);
p.addParameter('responseExtinctionMilliseconds', 20, @isnumeric);
p.addParameter('computeMosaic',true,@islogical);

% MOSAIC OPTIONS
p.addParameter('integrationTime', 5.0/1000, @isnumeric);

% HEX MOSAIC OPTIONS
p.addParameter('sConeMinDistanceFactor', 3.0, @isnumeric); % min distance between neighboring S-cones = f * local cone separation - to make the S-cone lattice semi-regular
p.addParameter('sConeFreeRadiusMicrons', 45, @isnumeric);
p.addParameter('latticeAdjustmentPositionalToleranceF', 0.01, @isnumeric);
p.addParameter('latticeAdjustmentDelaunayToleranceF', 0.001, @isnumeric);
p.addParameter('marginF', 1/sqrt(2.0), @isnumeric);     

% DIAGNOSTIC OPTIONS
p.addParameter('displayTrialBlockPartitionDiagnostics', true, @islogical);
p.addParameter('displayResponseComputationProgress', false, @islogical);

% RESPONSE MAP VISUALIZATION OPTIONS
p.addParameter('visualizeMosaic',true, @islogical); 
p.addParameter('visualizeOIsequence', true, @islogical);
p.addParameter('visualizeMosaicWithFirstEMpath', true, @islogical);
p.addParameter('visualizeSpatialScheme', true, @islogical);
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
p.addParameter('thresholdMethod', 'mlpt', @(x)ismember(x, {'svm', 'svmSpaceTimeSeparable', 'svmGaussianRF', 'mlpt', 'mlptGaussianRF','mlpe'}));
p.addParameter('thresholdSignal', 'isomerizations', @(x)ismember(x, {'isomerizations', 'photocurrents'}));
p.parse(paramsList);

userParams = p.Results;

end

% Update run params
function rParams = updateRunParams(rParams, userParams, backgroundLumIndex,stimDiamIndex)
% Update the spatial params
rParams = updateSpatialParams(rParams, userParams, stimDiamIndex); 

% Update the background wavelength
rParams.backgroundParams.backgroundWavelengthsNm = [userParams.wavelength];

% Scale radiance to produce desired background levels
rParams = updateBackgroundAndSpotParams(rParams, userParams, backgroundLumIndex);

% Update the temporal params
rParams = updateTemporalParams(rParams, userParams);

% Blur
rParams.oiParams = modifyStructParams(rParams.oiParams, ...
    'blur', userParams.blur, ...
    'pupilDiamMm', userParams.pupilDiamMm, ...
    'opticsModel', userParams.opticsModel);

% Update mosaic params
rParams = updateMosaicParams(rParams, userParams);
end


function rParams = updateMosaicParams(rParams, userParams)
%
rParams.mosaicParams = modifyStructParams(rParams.mosaicParams, ...
    'fieldOfViewDegs', userParams.backgroundSizeDegs, ... 
    'marginF', userParams.marginF, ...
    'innerSegmentSizeMicrons',userParams.innerSegmentSizeMicrons, ...
    'apertureBlur',userParams.apertureBlur, ...
    'mosaicRotationDegs',userParams.mosaicRotationDegs,...
    'coneDarkNoiseRate',userParams.coneDarkNoiseRate,...
    'LMSRatio',userParams.LMSRatio,...
    'coneSpacingMicrons', userParams.coneSpacingMicrons, ...
    'conePacking', userParams.conePacking, ...
    'integrationTimeInSeconds', userParams.integrationTime, ...
    'isomerizationNoise', 'random',...              % select from {'random', 'frozen', 'none'}
    'osNoise', 'random', ...                        % select from {'random', 'frozen', 'none'}
    'osModel', 'Linear');
        
% Additional params for eccentricity-based mosaics
if strcmp(userParams.conePacking, 'hex')
    rParams.mosaicParams = modifyStructParams(rParams.mosaicParams, ...
        'sConeMinDistanceFactor', userParams.sConeMinDistanceFactor, ...
        'sConeFreeRadiusMicrons', userParams.sConeFreeRadiusMicrons, ...
        'latticeAdjustmentPositionalToleranceF', userParams.latticeAdjustmentPositionalToleranceF, ...
        'latticeAdjustmentDelaunayToleranceF', userParams.latticeAdjustmentDelaunayToleranceF ...
        );
end
     
% Make sure mosaic noise parameters match the frozen noise flag
% passed in.  Doing this here means that things work on down the
% chain, while if we don't then there is a problem because routines
% create the wrong directory.
if (userParams.freezeNoise)
    if (strcmp(rParams.mosaicParams.isomerizationNoise, 'random'))
        rParams.mosaicParams.isomerizationNoise = 'frozen';
    end
    if (strcmp(rParams.mosaicParams.osNoise, 'random'))
        rParams.mosaicParams.osNoise = 'frozen';
    end
end
end


function rParams = updateTemporalParams(rParams, userParams)

stimulusDurationInSeconds = userParams.stimDurationMs/1000;
stimulusSamplingIntervalInSeconds = 1/userParams.stimRefreshRateHz;
windowTauInSeconds = nan; % square-wave
                
% Pre-stimulus time: allow photocurrent response to stabilize
responseStabilizationSeconds = ceil(userParams.responseStabilizationMilliseconds/1000/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;

% Post-stimulus time: allow photocurrent response to return to baseline
responseExtinctionSeconds = ceil(userParams.responseExtinctionMilliseconds/1000/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;

secondsToInclude = responseStabilizationSeconds+stimulusDurationInSeconds+responseExtinctionSeconds;
rParams.temporalParams = modifyStructParams(rParams.temporalParams, ...
    'frameRate', userParams.stimRefreshRateHz, ...
    'windowTauInSeconds', windowTauInSeconds, ...
    'stimulusSamplingIntervalInSeconds', stimulusSamplingIntervalInSeconds, ...
    'stimulusDurationInSeconds', stimulusDurationInSeconds, ...
    'secondsToInclude', secondsToInclude, ...
    'secondsForResponseStabilization', responseStabilizationSeconds, ...
    'secondsForResponseExtinction', responseExtinctionSeconds, ...
    'secondsToIncludeOffset', 0/1000, ...
    'emPathType', userParams.emPathType ...
);
end



function rParams = updateSpatialParams(rParams, userParams, stimDiamIndex)
rParams.spatialParams.spotSizeDegs       = userParams.spotDiametersMinutes(stimDiamIndex)/60;
rParams.spatialParams.backgroundSizeDegs = userParams.backgroundSizeDegs;
rParams.spatialParams.fieldOfViewDegs    = 1.1*userParams.backgroundSizeDegs;
rParams.spatialParams.row                = userParams.imagePixels;
rParams.spatialParams.col                = userParams.imagePixels;
rParams.spatialParams.viewingDistance    = 7.16;  
end

function rParams = updateBackgroundAndSpotParams(rParams, userParams, backgroundLumIndex)

%% Background params
deltaWl = 10;
S = [userParams.wavelength deltaWl 2];
wls = SToWls(S);
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
radiancePerUW = AOMonochromaticCornealPowerToRadiance(wls,rParams.backgroundParams.backgroundWavelengthsNm,1,rParams.oiParams.pupilDiamMm,rParams.spatialParams.backgroundSizeDegs^2);
xyzPerUW = T_xyz*radiancePerUW*(wls(2)-wls(1));
desiredLumCdM2 = userParams.luminances(backgroundLumIndex);
rParams.backgroundParams.backgroundCornealPowerUW = desiredLumCdM2/xyzPerUW(2);

%% Spot params
% Scale spot parameters into what we think a reasonable range is.
% In Davila & Geisler, they give thresholds in what they call
% threshold energy, which is spotLuminance*durationSecs*areaMinutes2.
%
% Just to get things into a reasonable scale, we'll find a
% luminance corresponding to a threshold energy of 10 for our
% smallest spot.

maxThresholdEnergy = 10;
minSpotAreaMin2 = pi*(min(userParams.spotDiametersMinutes)/2)^2;
rParams.maxSpotLuminanceCdM2 = maxThresholdEnergy/((userParams.stimDurationMs/1000)*minSpotAreaMin2);
rParams.colorModulationParams.startWl = userParams.wavelength;
rParams.colorModulationParams.endWl = userParams.wavelength+deltaWl;
rParams.colorModulationParams.deltaWl = deltaWl;
rParams.colorModulationParams.spotWavelengthNm = userParams.wavelength;
rParams.colorModulationParams.spotCornealPowerUW = rParams.maxSpotLuminanceCdM2/xyzPerUW(2);
rParams.colorModulationParams.contrast = 1;
end
    
