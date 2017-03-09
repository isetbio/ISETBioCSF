function [validationData, extraData, detectionThresholdContrast] = c_PoirsonAndWandell96Replicate(varargin)
%
% Compute thresholds to replicate the Poirson and Wandell 96 paper.
% We are finding thresholds for spatio-temporal stimuli defined in the full 3D LMS space.
%
%
% Key/value pairs
%    DATADIR OPTIONS
%       'useScratchTopLevelDirName'- true/false (default false). 
%           When true, the top level output directory is [scratch]. 
%           When false, it is the name of this script.
%
%    STIMULUS OPTIONS
%       'spatialFrequency': stimulus spatial frequency (default: 4 cpd)
%       'meanLuminance': stimulus mean luminance (default: 200 cd/m2)
%       'imagePixels' : how many pixels to use to represent the input stimulus
%
%    RESPONSE COMPUTATION OPTIONS
%       'nTrainingSamples' - how many response instances to compute. default: 128
%       'emPathType' - choose from {'none', 'frozen', 'random'}. Type of emPath: 
%           'none'   : results to 1 eye movement at (0,0)
%           'frozen' : results to the identical emPath being applied to each computed instance
%           'frozen0': all zeros path
%           'random' : results in a different emPath being applied to each computed instance
%       'parforWorkersNum' - 0 .. 12 (default: 12). How many workers to use for the computations.
%         use 0: for a serial for loop
%         use > 0: for a parfor loop with desired number of workers
%
%     MOSAIC OPTIONS
%       'coneMosaicPacking' - choose from {'hex' (default), 'hexReg', 'rect'}
%       'coneMosaicFOVDegs - the spatial extent of the cone mosaic (default: 1.00 degs)     
%       'freezeNoise' - true/false (default true).  Freezes isomerization and photocurrent noise so that results are reproducible
%       'computeResponses' - true/false (default true).  Do the computations.
%       'computeMosaic' - true/false (default true). Compute a cone mosaic or load one (good for large hex mosaics which take a while to compute)
%       'visualizeMosaic' - true/false (default true). Wether to visualize the cone mosaic
%
%     DIAGNOSTIC OPTIONS
%       'displayTrialBlockPartitionDiagnostics', true/false (default true). Wether to display trial block diagnostics.
%       'displayResponseComputationProgress', true/false. (default false) Wether to display the response computation progress
%
%     GENERAL VISUALIZATION OPTIONS
%       'generatePlots' - true/false (default false).  Produce response
%           visualizations.  Set to false when running big jobs on clusters or
%           in parfor loops, as plotting doesn't seem to play well with those conditions.
%
%     RESPONSE MAP VISUALIZATION OPTIONS
%       'visualizeResponses' - true/false (default true). Call the fancy visualize response routine.
%       'visualizeSpatialScheme' - true/false (default false). Visualize the relationship between mosaic and stimulus.
%       'visualizationFormat' - How to arrange visualized response maps.  Available options: 
%            'montage' (default). Create a montage of all frames 
%            'video'.             Animate through all frames. 
%       'visualizedResponseNormalization' - How to normalize visualized response maps. Available options: 
%           'submosaicBasedZscore', 
%           'LMSabsoluteResponseBased', 
%           'LMabsoluteResponseBased', 
%           'MabsoluteResponseBased'
%
%       'visualizePerformance' - true/false (default true). Plot the performance function
%
%     PERFORMANCE COMPUTATION OPTIONS
%       'performanceSignal' - either 'isomerizations', or 'photocurrents'
%       'performanceClassifier' - 'svm', 'mlpt', 'mlpe' (latter 2 applicable to Poisson noise, i.e. absorptions only)
%       'performanceClassifierTrainingSamples' - how many training samples
%       to use, default [], i.e. all
%       'performanceEvidenceIntegrationTime' - How many time bins should the classifier use
%       'pcaComponentsNum' - how many PCA components to use in the SVM methods - default: 60
%            (default: [], which corresponds to all the available time bins)
%       'findPerformance' - true/false (default true).  Find performance.
%       'fitPsychometric' - true/false (default true).  Fit psychometric functions.
%
%

%% Parse input
p = inputParser;
p.addParameter('useScratchTopLevelDirName', false, @islogical);
% STIMULUS OPTIONS
p.addParameter('spatialFrequency', 4.0, @isnumeric);
p.addParameter('meanLuminance', 200, @isnumeric);
p.addParameter('imagePixels', 256, @isnumeric);
p.addParameter('nContrastsPerDirection', 16, @isnumeric);
p.addParameter('lowContrast', 0.0001, @isnumeric);
p.addParameter('highContrast', 0.2, @isnumeric);
% RESPONSE COMPUTATION OPTIONS
p.addParameter('nTrainingSamples',128, @isnumeric);
p.addParameter('emPathType','frozen',@(x)ismember(x, {'none', 'frozen', 'frozen0', 'random'}));
p.addParameter('computeResponses',true,@islogical);
p.addParameter('computeMosaic',false,@islogical);
p.addParameter('parforWorkersNum', 12, @isnumeric);
% MOSAIC OPTIONS
p.addParameter('freezeNoise',true, @islogical);
p.addParameter('coneMosaicPacking', 'hex', @(x)ismember(x, {'hex', 'hexReg', 'rect'}));
p.addParameter('coneMosaicFOVDegs', 1.0, @isnumeric);
% DIAGNOSTIC OPTIONS
p.addParameter('displayTrialBlockPartitionDiagnostics', true, @islogical);
p.addParameter('displayResponseComputationProgress', false, @islogical);
% VISUALIZATION OPTIONS
p.addParameter('generatePlots',true,@islogical);
% RESPONSE MAP VISUALIZATION OPTIONS
p.addParameter('visualizeMosaic',true, @islogical); 
p.addParameter('visualizeResponses', false,@islogical);
p.addParameter('visualizeSpatialScheme', false, @islogical);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.addParameter('visualizationFormat', 'montage', @ischar);
p.addParameter('visualizePerformance', false, @islogical);
% PERFORMANCE COMPUTATION OPTIONS
p.addParameter('performanceClassifier', 'svm', @(x)ismember(x, {'svm', 'svmSpaceTimeSeparable', 'svmV1FilterBank', 'mlpt', 'mlpe'}));
p.addParameter('performanceSignal', 'isomerizations', @(x)ismember(x, {'isomerizations', 'photocurrents'}));
p.addParameter('performanceClassifierTrainingSamples', [], @isnumeric);
p.addParameter('performanceEvidenceIntegrationTime', [], @isnumeric);
p.addParameter('pcaComponentsNum', 60, @isnumeric);
p.addParameter('findPerformance',true,@islogical);
p.addParameter('fitPsychometric',true,@islogical);
p.parse(varargin{:});

% Ensure visualizationFormat has a valid value
visualizationFormat = p.Results.visualizationFormat;
if (strcmp(visualizationFormat, 'montage')) || (strcmp(visualizationFormat, 'video'))
else
    error('visualizationFormat must be set to either ''montage'' or ''video''. Current value: ''%s''.', visualizationFormat);
end

% Returned arguments
validationData = [];
extraData = [];
detectionThresholdContrast = [];

% Start with default
rParams = responseParamsGenerate;

%% Set the  topLevelDir name
if (~p.Results.useScratchTopLevelDirName)
    rParams.topLevelDirParams.name = mfilename;
end

% Modify spatial params to match P&W '96 - constant cycles condition
cyclesBandwidthProduct = 3.8;
gaussianFWHMDegs = cyclesBandwidthProduct/p.Results.spatialFrequency;
    
rParams.spatialParams = modifyStructParams(rParams.spatialParams, ...
        'windowType', 'Gaussian', ...
        'cyclesPerDegree', p.Results.spatialFrequency, ...
        'gaussianFWHMDegs', gaussianFWHMDegs, ...
        'fieldOfViewDegs',3*gaussianFWHMDegs, ...             % In P&W 1996, in the constant cycle condition, this was 10 deg (Section 2.2, p 517)
        'viewingDistance', 0.75, ...            % vd in meters
        'ang', 0,  ...                          % orientation in radians
        'ph', 0, ...                            % spatial phase in radians
        'row', p.Results.imagePixels, ...
        'col', p.Results.imagePixels);
  
% Modify background params to match P&W '96 (536.2 in paper)
luminancePW96 = p.Results.meanLuminance;  % limit to 200 for now because the photocurrent model is validated up to this luminance level
baseLum = 50;

rParams.backgroundParams = modifyStructParams(rParams.backgroundParams, ...
    'backgroundxyY', [0.38 0.39 baseLum]',...
    'monitorFile', 'CRT-MODEL', ...
    'leakageLum', 1.0, ...
    'lumFactor', luminancePW96/baseLum);

% Modify temporal params to match P&W'96
frameRate = 87;                                     % their CRT had 87 Hz refresh rate
windowTauInSeconds = 165/1000;
stimulusSamplingIntervalInSeconds = 1/frameRate;
stimulusDurationInSeconds = 4.5*windowTauInSeconds;
% Allow around 100 milliseconds for response to stabilize
responseStabilizationSeconds = ceil(100/1000/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
rParams.temporalParams = modifyStructParams(rParams.temporalParams, ...
    'frameRate', frameRate, ...
    'windowTauInSeconds', windowTauInSeconds, ...
    'stimulusSamplingIntervalInSeconds', stimulusSamplingIntervalInSeconds, ...
    'stimulusDurationInSeconds', stimulusDurationInSeconds, ...
    'secondsForResponseStabilization', 60/1000, ...
    'secondsToInclude', stimulusDurationInSeconds, ...
    'secondsForResponseStabilization', responseStabilizationSeconds, ...
    'secondsToIncludeOffset', 0/1000, ...
    'emPathType', p.Results.emPathType ...
);
   
equivalentIntegrationTime = [];
if strcmp(rParams.temporalParams.emPathType, 'none')
    % No eye movements
    % Adjust the integrationTime to match the area under the curve of the
    % stimulus modulation function
    [stimulusTimeAxis, stimulusModulationFunction, ~] = gaussianTemporalWindowCreate(rParams.temporalParams);
    timeIndicesToKeep = find(abs(stimulusTimeAxis-rParams.temporalParams.secondsToIncludeOffset) <= rParams.temporalParams.secondsToInclude/2);
    areaUnderTheCurve = sum(stimulusModulationFunction(timeIndicesToKeep)) * (stimulusTimeAxis(2)-stimulusTimeAxis(1));
    equivalentIntegrationTime = areaUnderTheCurve;
    
    % Run with a single exposure == windowTau and no eye movements
    % Equate stimulusSamplingIntervalInSeconds to stimulusDurationInSeconds to generate 1 time point only.
    rParams.temporalParams = modifyStructParams(rParams.temporalParams, ...
        'stimulusDurationInSeconds', stimulusDurationInSeconds, ...
        'stimulusSamplingIntervalInSeconds',  stimulusDurationInSeconds, ... 
        'secondsToInclude', stimulusDurationInSeconds ...  
    );
end

% Modify mosaic parameters
rParams.mosaicParams = modifyStructParams(rParams.mosaicParams, ...
    'conePacking', p.Results.coneMosaicPacking, ...                       
    'fieldOfViewDegs', p.Results.coneMosaicFOVDegs, ... 
    'integrationTimeInSeconds', 6/1000, ...
    'isomerizationNoise', 'random',...               % select from {'random', 'frozen', 'none'}
    'osNoise', 'random', ...                         % select from {'random', 'frozen', 'none'}
    'osModel', 'Linear');

if (p.Results.freezeNoise)
    if (strcmp(rParams.mosaicParams.isomerizationNoise, 'random'))
        rParams.mosaicParams.isomerizationNoise = 'frozen';
    end
    if (strcmp(rParams.mosaicParams.osNoise, 'random'))
        rParams.mosaicParams.osNoise = 'frozen';
    end
end

% Parameters that define the LMS instances we'll generate
% Here, we are generating an L+M stimulus (azimuth = 45, elevation = 0);
testDirectionParams = instanceParamsGenerate('instanceType', 'LMSPlane');
testDirectionParams = modifyStructParams(testDirectionParams, ...
    'trialsNum', p.Results.nTrainingSamples, ...
    'startAzimuthAngle', 45, ...
    'nAzimuthAngles', 1, ...
    'startElevationAngle', 0, ...
    'nElevationAngles', 1, ...
    'nContrastsPerDirection', p.Results.nContrastsPerDirection, ...
    'lowContrast', p.Results.lowContrast, ...
    'highContrast', p.Results.highContrast, ...
    'contrastScale', 'log' ...    % choose between 'linear' and 'log'  
);

%% Compute response instances
if (p.Results.computeResponses)
    tBegin = clock;
    t_coneCurrentEyeMovementsResponseInstances(...
          'rParams',rParams,...
          'testDirectionParams',testDirectionParams,...
          'centeredEMPaths',true, ...
          'freezeNoise', p.Results.freezeNoise, ...
          'compute',p.Results.computeResponses, ...
          'computeMosaic', p.Results.computeMosaic, ... 
          'visualizeMosaic', p.Results.visualizeMosaic, ...
          'parforWorkersNum', p.Results.parforWorkersNum, ...  % no more than these many workers
          'overrideMosaicIntegrationTime', equivalentIntegrationTime, ... 
          'displayTrialBlockPartitionDiagnostics', p.Results.displayTrialBlockPartitionDiagnostics, ...
          'generatePlots', p.Results.generatePlots, ...
          'visualizeSpatialScheme', p.Results.visualizeSpatialScheme, ...
          'visualizeResponses', p.Results.visualizeResponses, ... 
          'visualizedResponseNormalization', p.Results.visualizedResponseNormalization, ...
          'visualizationFormat', p.Results.visualizationFormat, ...
          'workerID', find(p.Results.displayResponseComputationProgress));
    tEnd = clock;
    timeLapsed = etime(tEnd,tBegin);
    fprintf('Computation of isomerization & photocurrent responses was completed in %f minutes. \n', timeLapsed/60);
end

%% Visualize response instances
if (p.Results.visualizeResponses) 
    t_coneCurrentEyeMovementsResponseInstances(...
          'rParams',rParams,...
          'testDirectionParams',testDirectionParams,...
          'centeredEMPaths',true, ...
          'freezeNoise', p.Results.freezeNoise, ...
          'compute', false, ...
          'computeMosaic', false, ... 
          'visualizedResponseNormalization', p.Results.visualizedResponseNormalization, ...
          'visualizationFormat', p.Results.visualizationFormat, ...
          'generatePlots', true, ...
          'visualizeResponses', true);
end % visualizeResponses


%% Find performance
% Parameters related to how we find thresholds from responses
% Use default
thresholdParams = thresholdParamsGenerate;
   
% Reduce # of trials used if computation is not feasible (ie. PCA)
% Here we are using all of them
if (isempty(p.Results.performanceClassifierTrainingSamples))
    thresholdParams.trialsUsed = p.Results.nTrainingSamples;
else
    thresholdParams.trialsUsed = p.Results.performanceClassifierTrainingSamples; % p.Results.nTrainingSamples;
end


thresholdParams = modifyStructParams(thresholdParams, ...
    'method', p.Results.performanceClassifier, ...
    'STANDARDIZE', false, ...
    'standardizeSVMpredictors', false, ...
    'PCAComponents', p.Results.pcaComponentsNum, ...
    'evidenceIntegrationTime', p.Results.performanceEvidenceIntegrationTime, ...
    'signalSource', p.Results.performanceSignal);

if (p.Results.findPerformance) || (p.Results.visualizePerformance)
    rParams.plotParams = modifyStructParams(rParams.plotParams, ...
        'axisFontSize', 12, ...
        'labelFontSize', 14, ...
        'lineWidth', 1.5);
    
    t_colorDetectFindPerformance(...
        'rParams',rParams, ...
        'testDirectionParams',testDirectionParams,...
        'thresholdParams',thresholdParams, ...
        'compute',p.Results.findPerformance, ...
        'visualizeSpatialScheme', p.Results.visualizeSpatialScheme, ...
        'plotSvmBoundary',false, ...
        'plotPsychometric',true ...
        );
    
    % Fit psychometric functions
    if (p.Results.fitPsychometric)
      d = t_plotDetectThresholdsOnLMPlane(...
          'rParams',rParams, ...
          'instanceParams',testDirectionParams, ...
          'thresholdParams',thresholdParams, ...
          'plotPsychometric',p.Results.visualizePerformance, ...
          'plotEllipse',false);
      
      detectionThresholdContrast = d.thresholdContrasts;
    end
        
end


        
end

