function [detectionThresholdContrast, pedestalLuminanceList] = c_TVIcurve(varargin)

%% Parse input
p = inputParser;
p.addParameter('useScratchTopLevelDirName', false, @islogical);
% STIMULUS OPTIONS
p.addParameter('imagePixels', 128, @isnumeric);
p.addParameter('fieldOfViewDegs', 0.8, @isnumeric);
p.addParameter('testDiameterDegs', 0.4, @isnumeric);
p.addParameter('nContrastsPerPedestalLuminance', 12, @isnumeric);
p.addParameter('lowContrast', 1e-4, @isnumeric);  % 3e-4
p.addParameter('highContrast', 0.1, @isnumeric);  % 0.3
p.addParameter('nPedestalLuminanceLevels', 10, @isnumeric);
p.addParameter('lowPedestalLuminance', 3.3, @isnumeric); 
p.addParameter('highPedestalLuminance', 125, @isnumeric);
p.addParameter('pedestalLuminanceListIndicesUsed', [], @isnumeric);
% RESPONSE COMPUTATION OPTIONS
p.addParameter('nTrainingSamples',512, @isnumeric);
p.addParameter('emPathType','frozen0',@(x)ismember(x, {'none', 'frozen', 'frozen0', 'random'}));
p.addParameter('computeResponses',true,@islogical);
p.addParameter('computeMosaic',false,@islogical);
p.addParameter('ramPercentageEmployed', 1.0, @isnumeric);
p.addParameter('parforWorkersNum', 12, @isnumeric);
% MOSAIC OPTIONS
p.addParameter('freezeNoise',true, @islogical);
p.addParameter('osNoise', 'random', @(x)ismember(x, {'random', 'frozen', 'none'}));
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
p.addParameter('visualizeOuterSegmentFilters', false,@islogical);
p.addParameter('visualizeSpatialScheme', false, @islogical);
p.addParameter('visualizeTransformedSignals',false, @islogical);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.addParameter('visualizationFormat', 'montage', @ischar);
p.addParameter('visualizePerformance', false, @islogical);
% PERFORMANCE COMPUTATION OPTIONS
p.addParameter('performanceClassifier', 'mlpt', @(x)ismember(x, {'svm', 'svmSpaceTimeSeparable', 'svmGaussianRF', 'mlpt', 'mlpe'}));
p.addParameter('performanceSignal', 'isomerizations', @(x)ismember(x, {'isomerizations', 'photocurrents'}));
p.addParameter('performanceClassifierTrainingSamples', [], @isnumeric);
p.addParameter('performanceEvidenceIntegrationTime', [], @isnumeric);
p.addParameter('pcaComponentsNum', 60, @isnumeric);
p.addParameter('findPerformance',true,@islogical);
p.addParameter('fitPsychometric',true,@islogical);
p.parse(varargin{:});


% Return params
detectionThresholdContrast = {};

% Ensure visualizationFormat has a valid value
visualizationFormat = p.Results.visualizationFormat;
if (strcmp(visualizationFormat, 'montage')) || (strcmp(visualizationFormat, 'video'))
else
    error('visualizationFormat must be set to either ''montage'' or ''video''. Current value: ''%s''.', visualizationFormat);
end


% Start with default params
rParams = responseParamsGenerate;

%% Set the  topLevelDir name
if (~p.Results.useScratchTopLevelDirName)
    rParams.topLevelDirParams.name = mfilename;
end

% Modify spatial params
rParams.spatialParams = spatialParamsGenerate('spatialType','pedestalDisk');
rParams.spatialParams = modifyStructParams(rParams.spatialParams, ...
        'testDiameterDegs', p.Results.testDiameterDegs, ...
        'fieldOfViewDegs', p.Results.fieldOfViewDegs, ... 
        'viewingDistance', 0.75, ...            % vd in meters
        'row', p.Results.imagePixels, ...
        'col', p.Results.imagePixels);
    

% Modify temporal params 
frameRate = 60;          
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
    'secondsToInclude', stimulusDurationInSeconds, ...
    'secondsForResponseStabilization', responseStabilizationSeconds, ...
    'secondsToIncludeOffset', 0/1000, ...
    'emPathType', p.Results.emPathType ...
);

% Modify mosaic parameters
rParams.mosaicParams = modifyStructParams(rParams.mosaicParams, ...
    'conePacking', p.Results.coneMosaicPacking, ...                       
    'fieldOfViewDegs', p.Results.coneMosaicFOVDegs, ... 
    'integrationTimeInSeconds', 6/1000, ...
    'isomerizationNoise', 'random',...               % select from {'random', 'frozen', 'none'}
    'osNoise', p.Results.osNoise, ...                % select from {'random', 'frozen', 'none'}
    'osModel', 'Linear');

if (p.Results.freezeNoise)
    if (strcmp(rParams.mosaicParams.isomerizationNoise, 'random'))
        rParams.mosaicParams.isomerizationNoise = 'frozen';
    end
    if (strcmp(rParams.mosaicParams.osNoise, 'random'))
        rParams.mosaicParams.osNoise = 'frozen';
    end
end

% Parameters that define the pedestal increments that we'll generate
testLuminanceParams = instanceParamsGenerate('instanceType', 'pedestalIncrements');
testLuminanceParams = modifyStructParams(testLuminanceParams, ...
    'trialsNum', p.Results.nTrainingSamples, ...
    'nContrastsPerDirection', p.Results.nContrastsPerPedestalLuminance, ...
    'lowContrast', p.Results.lowContrast, ...
    'highContrast', p.Results.highContrast, ...
    'contrastScale', 'log' ...    % choose between 'linear' and 'log'  
);


% Parameters that determine how we measure detection performance
% Start with default
thresholdParams = thresholdParamsGenerate;
   
% Reduce # of trials used if computation is not feasible (ie. PCA)
if (isempty(p.Results.performanceClassifierTrainingSamples))
    thresholdParams.trialsUsed = p.Results.nTrainingSamples;
else
    thresholdParams.trialsUsed = p.Results.performanceClassifierTrainingSamples;
end


thresholdParams = modifyStructParams(thresholdParams, ...
    'method', p.Results.performanceClassifier, ...
    'STANDARDIZE', false, ...
    'standardizeSVMpredictors', false, ...
    'PCAComponents', p.Results.pcaComponentsNum, ...
    'evidenceIntegrationTime', p.Results.performanceEvidenceIntegrationTime, ...
    'signalSource', p.Results.performanceSignal);


% Start timing
tBegin = clock;

% Loop over pedestal luminances
pedestalLuminanceList = logspace(log10(p.Results.lowPedestalLuminance), log10(p.Results.highPedestalLuminance), p.Results.nPedestalLuminanceLevels);
if (~isempty(p.Results.pedestalLuminanceListIndicesUsed))
    pedestalLuminanceList = pedestalLuminanceList(p.Results.pedestalLuminanceListIndicesUsed);
end
    
for pedestalLuminanceIndex = 1:numel(pedestalLuminanceList)
    
    %% Background params
    pedestalLuminance = pedestalLuminanceList(pedestalLuminanceIndex);
    baseLum = 50;
    rParams.backgroundParams = modifyStructParams(rParams.backgroundParams, ...
        'backgroundxyY', [0.31 0.316 baseLum]',...  % Illuminant-C chromaticity
        'monitorFile', 'CRT-MODEL', ...
        'leakageLum', 1.0, ...
        'lumFactor', pedestalLuminance/baseLum);

    %% Compute response instances
    if ((p.Results.computeResponses) || (p.Results.visualizeSpatialScheme)) 
        t_coneCurrentEyeMovementsResponseInstances(...
              'rParams',rParams,...
              'testDirectionParams',testLuminanceParams,...
              'centeredEMPaths',true, ...
              'freezeNoise', p.Results.freezeNoise, ...
              'compute',p.Results.computeResponses, ...
              'computeMosaic', p.Results.computeMosaic, ... 
              'visualizeMosaic', p.Results.visualizeMosaic, ...
              'ramPercentageEmployed', p.Results.ramPercentageEmployed, ...
              'parforWorkersNum', p.Results.parforWorkersNum, ...  % no more than these many workers
              'displayTrialBlockPartitionDiagnostics', p.Results.displayTrialBlockPartitionDiagnostics, ...
              'generatePlots', p.Results.generatePlots, ...
              'visualizeSpatialScheme', p.Results.visualizeSpatialScheme, ...
              'visualizeResponses', p.Results.visualizeResponses, ... 
              'visualizedResponseNormalization', p.Results.visualizedResponseNormalization, ...
              'visualizationFormat', p.Results.visualizationFormat, ...
              'workerID', find(p.Results.displayResponseComputationProgress));
    end % if (p.Results.computeResponses)
    
    %% Visualize response instances
    if ((p.Results.visualizeResponses) || (p.Results.visualizeOuterSegmentFilters))
    [~,~, osImpulseResponses] = t_coneCurrentEyeMovementsResponseInstances(...
          'rParams',rParams,...
          'testDirectionParams',testLuminanceParams,...
          'freezeNoise', p.Results.freezeNoise, ...
          'compute', false, ...
          'computeMosaic', false, ... 
          'visualizedResponseNormalization', p.Results.visualizedResponseNormalization, ...
          'visualizationFormat', p.Results.visualizationFormat, ...
          'generatePlots', true, ...
          'visualizeResponses', p.Results.visualizeResponses, ...
          'visualizeSpatialScheme', p.Results.visualizeSpatialScheme, ...
          'visualizeOuterSegmentFilters', p.Results.visualizeOuterSegmentFilters);
      
      if (p.Results.visualizeOuterSegmentFilters)
            irAmplitudes(pedestalLuminanceIndex,:) = squeeze(max(osImpulseResponses, [], 1));
      end
    end % visualizeResponses

    %% Compute/Visualize performance
    if (p.Results.findPerformance) || (p.Results.visualizePerformance) || (p.Results.fitPsychometric) || (p.Results.visualizeSpatialScheme)
        rParams.plotParams = modifyStructParams(rParams.plotParams, ...
            'axisFontSize', 12, ...
            'labelFontSize', 14, ...
            'lineWidth', 1.5);

        if (p.Results.findPerformance) || (p.Results.visualizePerformance) || (p.Results.visualizeSpatialScheme)
            t_colorDetectFindPerformance(...
                'rParams',rParams, ...
                'testDirectionParams',testLuminanceParams,...
                'thresholdParams',thresholdParams, ...
                'compute',p.Results.findPerformance, ...
                'visualizeSpatialScheme', p.Results.visualizeSpatialScheme, ...
                'visualizeTransformedSignals', p.Results.visualizeTransformedSignals, ...
                'plotSvmBoundary',false, ...
                'plotPsychometric',true ...
                );
        end
        
        % Fit psychometric functions
        if (p.Results.fitPsychometric) && ((p.Results.findPerformance) || (p.Results.visualizePerformance))
          d = t_plotDetectThresholds(...
              'rParams',rParams, ...
              'instanceParams',testLuminanceParams, ...
              'thresholdParams',thresholdParams, ...
              'plotPsychometric',p.Results.visualizePerformance);
          detectionThresholdContrast{pedestalLuminanceIndex} = d.thresholdContrasts;
        end
    end % if (p.Results.findPerformance) || (p.Results.visualizePerformance)
end % pedestalLuminanceIndex

tEnd = clock;
timeLapsed = etime(tEnd,tBegin);
fprintf('Full computation was completed in %f minutes. \n', timeLapsed/60, pedestalLuminance);

if (p.Results.visualizeOuterSegmentFilters)
    hFig = figure(999); clf;
    set(hFig, 'Position', [10 10 450 450], 'Color', [1 1 1]);
    plot(log10(pedestalLuminanceList), log10(irAmplitudes(:,1)), 'rs-', 'LineWidth', 1.5, 'MarkerSize', 14, 'MarkerFaceColor', [0.6 0.6 0.6]);
    hold on
    plot(log10(pedestalLuminanceList), log10(irAmplitudes(:,2)), 'gs-', 'LineWidth', 1.5, 'MarkerSize', 14, 'MarkerFaceColor', [0.6 0.6 0.6]);
    plot(log10(pedestalLuminanceList), log10(irAmplitudes(:,3)), 'bs-', 'LineWidth', 1.5, 'MarkerSize', 14, 'MarkerFaceColor', [0.6 0.6 0.6]);
    set(gca, 'FontSize', 14);
    set(gca, 'XLim', [log10(pedestalLuminanceList(1)*0.95) log10(pedestalLuminanceList(end)*1.05)]);
    grid on;
    box off;
    axis 'square'
    xlabel('log pedestal luminance (cd/m2)', 'FontWeight', 'bold');
    ylabel('log outer segment impulse response amplitude',  'FontWeight', 'bold');
    
    paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams};
    theProgram = mfilename;
    rwObject = IBIOColorDetectReadWriteBasic;
    data = 0;
    rwObject.write('os_ImpulseResponseSummary', data, paramsList, theProgram, ...
        'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
    
end


if (p.Results.findPerformance) || (p.Results.visualizePerformance)
    detectionThresholdContrast = cell2mat(detectionThresholdContrast);

    hFig = figure(1000); clf;
    set(hFig, 'Position', [10 10 450 900], 'Color', [1 1 1]);
    
    subplot(2,1,1)
    plot(log10(pedestalLuminanceList), log10(detectionThresholdContrast), 'rs-', 'LineWidth', 1.5, 'MarkerSize', 14, 'MarkerFaceColor', [1 0.8 0.8]);
    set(gca, 'XLim', [log10(pedestalLuminanceList(1)*0.95) log10(pedestalLuminanceList(end)*1.05)], 'YLim', [-3.4 -0.9]);
    xlabel('log pedestal luminance (cd/m2)', 'FontWeight', 'bold');
    ylabel('log threshold contrast', 'FontWeight', 'bold');
    title(sprintf('threshold contrast for \n%s (%s)\n',thresholdParams.signalSource, thresholdParams.method));
    set(gca, 'FontSize', 14);
    grid on;
    box off;
    axis 'square'
    
    
    subplot(2,1,2);
    thresholdDeltaLuminance = pedestalLuminanceList .* detectionThresholdContrast;
    plot(log10(pedestalLuminanceList), log10(thresholdDeltaLuminance), 'rs-', 'LineWidth', 1.5, 'MarkerSize', 14, 'MarkerFaceColor', [1 0.8 0.8]);
    set(gca, 'XLim', [log10(pedestalLuminanceList(1)*0.95) log10(pedestalLuminanceList(end)*1.05)], 'YLim', [-2.0 0.2]);
    title(sprintf('threshold delta-luminance for \n%s (%s)\n',thresholdParams.signalSource, thresholdParams.method))
    xlabel('log pedestal luminance (cd/m2)', 'FontWeight', 'bold');
    ylabel('log threshold delta luminance', 'FontWeight', 'bold');
    set(gca, 'FontSize', 14);
    grid on;
    box off;
    axis 'square'
    drawnow

    
    paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams, thresholdParams};
    theProgram = mfilename;
    rwObject = IBIOColorDetectReadWriteBasic;
    data = 0;
    rwObject.write('c_TCIcurveSummary', data, paramsList, theProgram, ...
        'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
end
       
end

