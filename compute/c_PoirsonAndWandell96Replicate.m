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
%     % RESPONSE COMPUTATION OPTIONS
%       'imagePixels' : how many pixels to use to represent the input stimulus
%       'nTrainingSamples' - how many response instances to compute. default: 128
%       'emPathType' - choose from {'none', 'frozen', 'random'}. Type of emPath: 
%           'none'   : results to 1 eye movement at (0,0)
%           'frozen' : results to the identical emPath being applied to each computed instance
%           'random' : results in a different emPath being applied to each computed instance
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
%       'visualizationFormat' - How to arrange visualized response maps.  Available options: 
%            'montage' (default). Create a montage of all frames 
%            'video'.             Animate through all frames. 
%       'visualizedResponseNormalization' - How to normalize visualized response maps. Available options: 
%           'submosaicBasedZscore', 
%           'LMSabsoluteResponseBased', 
%           'LMabsoluteResponseBased', 
%           'MabsoluteResponseBased'
%
%     PERFORMANCE COMPUTATION OPTIONS
%       'performanceSignal' - either 'isomerizations', or 'photocurrents'
%       'findPerformance' - true/false (default true).  Find performance.
%       'fitPsychometric' - true/false (default true).  Fit psychometric functions.
%
%

%% Parse input
p = inputParser;
p.addParameter('useScratchTopLevelDirName', false, @islogical);
% RESPONSE COMPUTATION OPTIONS
p.addParameter('imagePixels',500, @isnumeric);
p.addParameter('nTrainingSamples',128, @isnumeric);
p.addParameter('emPathType','frozen',@(x)ismember(x, {'none', 'frozen', 'random'}));
p.addParameter('freezeNoise',true,@islogical);
p.addParameter('computeResponses',true,@islogical);
p.addParameter('computeMosaic',false,@islogical);
% DIAGNOSTIC OPTIONS
p.addParameter('displayTrialBlockPartitionDiagnostics', true, @islogical);
% VISUALIZATION OPTIONS
p.addParameter('generatePlots',true,@islogical);
% RESPONSE MAP VISUALIZATION OPTIONS
p.addParameter('visualizeResponses',true,@islogical);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.addParameter('visualizationFormat', 'montage', @ischar);
% PERFORMANCE COMPUTATION OPTIONS
p.addParameter('performanceSignal', 'isomerizations', @ischar);
p.addParameter('findPerformance',true,@islogical);
p.addParameter('fitPsychometric',true,@islogical);

p.parse(varargin{:});

% Ensure visualizationFormat has a valid value
visualizationFormat = p.Results.visualizationFormat;
if (strcmp(visualizationFormat, 'montage')) || (strcmp(visualizationFormat, 'video'))
else
    error('visualizationFormat must be set to either ''montage'' or ''video''. Current value: ''%s''.', visualizationFormat);
end

% Ensure the performanceSignal has a valid value
performanceSignal = p.Results.performanceSignal;
if (strcmp(visualizationFormat, 'isomerizations')) || (strcmp(visualizationFormat, 'photocurrents'))
else
    error('performanceSignal must be set to either ''isomerizations'' or ''photocurrents''. Current value: ''%s''.', performanceSignal);
end


% Start with default
rParams = responseParamsGenerate;

%% Set the  topLevelDir name
if (~p.Results.useScratchTopLevelDirName)
    rParams.topLevelDirParams.name = mfilename;
end

% Modify spatial params to match P&W '96
rParams.spatialParams = modifyStructParams(rParams.spatialParams, ...
        'windowType', 'Gaussian', ...
        'cyclesPerDegree', 4.0, ...
        'gaussianFWHMDegs', 1.9, ...
        'fieldOfViewDegs', 10.0, ...             % In P&W 1996, in the constant cycle condition, this was 10 deg (Section 2.2, p 517)
        'viewingDistance', 0.75, ...            % vd in meters
        'ang', 0,  ...                          % orientation in radians
        'ph', 0, ...                            % spatial phase in radians
        'row', p.Results.imagePixels, ...
        'col', p.Results.imagePixels);
  
% Modify background params to match P&W '96
luminancePW96 = 536.2;
luminancePW96 = 200.0;  % limit to 200 for now because the photocurrent model is validated up to this luminance level
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
stimulusDurationInSeconds = 1.5*windowTauInSeconds;
rParams.temporalParams = modifyStructParams(rParams.temporalParams, ...
    'frameRate', frameRate, ...
    'windowTauInSeconds', windowTauInSeconds, ...
    'stimulusSamplingIntervalInSeconds', stimulusSamplingIntervalInSeconds, ...
    'stimulusDurationInSeconds', stimulusDurationInSeconds, ...
    'secondsToInclude', 300/1000, ...
    'secondsToIncludeOffset', 0/1000, ...
    'emPathType', p.Results.emPathType ...
);
    
if strcmp(p.Results.emPathType, 'none')
    % No eye movements
    [stimulusTimeAxis, stimulusModulationFunction, ~] = gaussianTemporalWindowCreate(temporalParams);
    timeIndicesToKeep = find(abs(stimulusTimeAxis-temporalParams.secondsToIncludeOffset) <= temporalParams.secondsToInclude/2);
    areaUnderTheCurve = sum(stimulusModulationFunction(timeIndicesToKeep)) * (stimulusTimeAxis(2)-stimulusTimeAxis(1));
    equivalentIntegrationTime = areaUnderTheCurve
    
    % Run with a single exposure == windowTau and no eye movements
    % Equate stimulusSamplingIntervalInSeconds to stimulusDurationInSeconds to generate 1 time point only.
    rParams.temporalParams = modifyStructParams(rParams.temporalParams, ...
        'stimulusDurationInSeconds', stimulusDurationInSeconds, ...
        'stimulusSamplingIntervalInSeconds',  stimulusDurationInSeconds, ... 
        'secondsToInclude', stimulusDurationInSeconds ...  
    );

    rParams.mosaicParams.integrationTime = equivalentIntegrationTime;
end

% Modify mosaic parameters
conePacking = 'hex';        % Hexagonal mosaic
%conePacking = 'rect';       % Rectangular mosaic
rParams.mosaicParams = modifyStructParams(rParams.mosaicParams, ...
    'conePacking', conePacking, ...                       
    'fieldOfViewDegs', rParams.spatialParams.fieldOfViewDegs*0.125, ... 
    'integrationTimeInSeconds', 6/1000, ...
    'isomerizationNoise', 'frozen',...              % select from {'random', 'frozen', 'none'}
    'osNoise', 'frozen', ...                        % select from {'random', 'frozen', 'none'}
    'osModel', 'Linear');
        
% Parameters that define the LMS instances we'll generate
% Here, we are generating an L+M stimulus (azimuth = 45, elevation = 0);
testDirectionParams = instanceParamsGenerate('instanceType', 'LMSPlane');
testDirectionParams = modifyStructParams(testDirectionParams, ...
    'trialsNum', p.Results.nTrainingSamples, ...
    'startAzimuthAngle', 45, ...
    'nAzimuthAngles', 1, ...
    'startElevationAngle', 0, ...
    'nElevationAngles', 1, ...
    'nContrastsPerDirection', 12, ...
    'lowContrast', 0.0001, ...
    'highContrast', 0.1, ...
    'contrastScale', 'log' ...    % choose between 'linear' and 'log'  
);

% Parameters related to how we find thresholds from responses
% Use default
thresholdParams = thresholdParamsGenerate;
   

%% Compute response instances
if (p.Results.computeResponses)
    tBegin = clock;
    t_coneCurrentEyeMovementsResponseInstances(...
          'rParams',rParams,...
          'testDirectionParams',testDirectionParams,...
          'compute',p.Results.computeResponses, ...
          'computeMosaic', p.Results.computeMosaic, ... 
          'freezeNoise', p.Results.freezeNoise, ...
          'trialBlocks', -1, ...                    % automatically decide trialBlocks based on system resources
          'displayTrialBlockPartitionDiagnostics', p.Results.displayTrialBlockPartitionDiagnostics, ...
          'generatePlots', p.Results.generatePlots, ...
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
          'compute', false, ...
          'computeMosaic', false, ... 
          'visualizedResponseNormalization', p.Results.visualizedResponseNormalization, ...
          'visualizationFormat', p.Results.visualizationFormat, ...
          'generatePlots', true, ...
          'visualizeResponses', true);
end % visualizeResponses


%% Find performance, template max likeli
thresholdParams.method = 'mlpt';
% Reduce # of trials used to make computation feasible
thresholdParams.trialsUsed = 128;
thresholdParams.signalSource = performanceSignal;

if (p.Results.findPerformance)
    t_colorDetectFindPerformance(...
        'rParams',rParams, ...
        'testDirectionParams',testDirectionParams,...
        'thresholdParams',thresholdParams, ...
        'compute',true, ...
        'plotSvmBoundary',false, ...
        'plotPsychometric',true ...
        );
end
        
end

