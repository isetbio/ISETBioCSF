function [validationData, extraData] = c_PoirsonAndWandell96Replicate(varargin)
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
%
%     MOSAIC OPTIONS
%       'coneMosaicPacking' - choose from {'hex' (default), 'hexReg', 'rect'}
%       'coneMosaicFOVDegs - the spatial extent of the cone mosaic (default: 1.00 degs)     
%       'freezeNoise' - true/false (default true).  Freezes isomerization and photocurrent noise so that results are reproducible
%       'computeResponses' - true/false (default true).  Do the computations.
%       'computeMosaic' - true/false (default true). Compute a cone mosaic or load one (good for large hex mosaics which take a while to compute)
%
%     DIAGNOSTIC OPTIONS
%       'displayTrialBlockPartitionDiagnostics', true/false. Wether to display trial block diagnostics.
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
%       'performanceEvidenceIntegrationTime' - How many time bins should the classifier use
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
% MOSAIC OPTIONS
p.addParameter('freezeNoise',true, @islogical);
p.addParameter('coneMosaicPacking', 'hex', @(x)ismember(x, {'hex', 'hexReg', 'rect'}));
p.addParameter('coneMosaicFOVDegs', 1.0, @isnumeric);
% DIAGNOSTIC OPTIONS
p.addParameter('displayTrialBlockPartitionDiagnostics', true, @islogical);
% VISUALIZATION OPTIONS
p.addParameter('generatePlots',true,@islogical);
% RESPONSE MAP VISUALIZATION OPTIONS
p.addParameter('visualizeResponses', false,@islogical);
p.addParameter('visualizeSpatialScheme', false, @islogical);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.addParameter('visualizationFormat', 'montage', @ischar);
p.addParameter('visualizePerformance', false, @islogical);
% PERFORMANCE COMPUTATION OPTIONS
p.addParameter('performanceClassifier', 'svm', @(x)ismember(x, {'svm', 'svmV1FilterBank', 'mlpt', 'mlpe'}));
p.addParameter('performanceSignal', 'isomerizations', @(x)ismember(x, {'isomerizations', 'photocurrents'}));
p.addParameter('performanceEvidenceIntegrationTime', [], @isnumeric);
p.addParameter('findPerformance',true,@islogical);
p.addParameter('fitPsychometric',true,@islogical);
p.parse(varargin{:});

% Ensure visualizationFormat has a valid value
visualizationFormat = p.Results.visualizationFormat;
if (strcmp(visualizationFormat, 'montage')) || (strcmp(visualizationFormat, 'video'))
else
    error('visualizationFormat must be set to either ''montage'' or ''video''. Current value: ''%s''.', visualizationFormat);
end

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
stimulusDurationInSeconds = 3.0*windowTauInSeconds;
% Allow around 100 milliseconds for response to stabilize
responseStabilizationSeconds = ceil(100/1000/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
rParams.temporalParams = modifyStructParams(rParams.temporalParams, ...
    'frameRate', frameRate, ...
    'windowTauInSeconds', windowTauInSeconds, ...
    'stimulusSamplingIntervalInSeconds', stimulusSamplingIntervalInSeconds, ...
    'stimulusDurationInSeconds', stimulusDurationInSeconds, ...
    'secondsForResponseStabilization', 60/1000, ...
    'secondsToInclude', 800/1000, ...
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
          'compute',p.Results.computeResponses, ...
          'computeMosaic', p.Results.computeMosaic, ... 
          'parforWorkersNum', 12, ...  % no more than 12 workers
          'overrideMosaicIntegrationTime', equivalentIntegrationTime, ... 
          'freezeNoise', p.Results.freezeNoise, ...
          'trialBlocks', -1, ...                    % automatically decide trialBlocks based on system resources
          'displayTrialBlockPartitionDiagnostics', p.Results.displayTrialBlockPartitionDiagnostics, ...
          'generatePlots', p.Results.generatePlots, ...
          'visualizeSpatialScheme', p.Results.visualizeSpatialScheme, ...
          'visualizeResponses', p.Results.visualizeResponses, ... 
          'visualizedResponseNormalization', p.Results.visualizedResponseNormalization, ...
          'visualizationFormat', p.Results.visualizationFormat, ...
          'workerID', 1);
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
thresholdParams.trialsUsed = p.Results.nTrainingSamples;

thresholdParams = modifyStructParams(thresholdParams, ...
    'method', p.Results.performanceClassifier, ...
    'STANDARDIZE', false, ...
    'standardizeSVMpredictors', false, ...
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
        'parforWorkersNum', 0, ... % do a serial for loop
        'plotSvmBoundary',false, ...
        'plotPsychometric',true ...
        );
end

end

