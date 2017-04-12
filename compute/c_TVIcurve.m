function [detectionThresholdContrast, pedestalLuminanceList] = c_TVIcurve(varargin)

%% Parse input
p = inputParser;
p.addParameter('useScratchTopLevelDirName', false, @islogical);
% STIMULUS OPTIONS - SPATIAL
p.addParameter('imagePixels', 128, @isnumeric);
p.addParameter('fieldOfViewDegs', 0.8, @isnumeric);
p.addParameter('testDiameterDegs', 0.4, @isnumeric);
% STIMULUS OPTIONS - TEMPORAL
p.addParameter('stimulusTemporalEnvelope', 'Gaussian', @(x)ismember(x, {'Gaussian', 'square'}));
p.addParameter('stimulusDurationSecs', 495/1000, @isnumeric);
% STIMULUS OPTIONS - CONTRAST/LUMINANCE
p.addParameter('nContrastsPerPedestalLuminance', 12, @isnumeric);
p.addParameter('lowContrast', 3e-3, @isnumeric);  
p.addParameter('highContrast', 0.1, @isnumeric); 
p.addParameter('nPedestalLuminanceLevels', 10, @isnumeric);
p.addParameter('lowPedestalLuminance', 3.3, @isnumeric);         % Fred's model developed in luminance range: 500-20,000 R*/cone/sec], which corresponds to luminance range: 3.3 - 125 cd/m2
p.addParameter('highPedestalLuminance', 125, @isnumeric);
p.addParameter('pedestalLuminanceListIndicesUsed', [], @isnumeric);
% RESPONSE COMPUTATION OPTIONS
p.addParameter('nTrainingSamples',512, @isnumeric);
p.addParameter('emPathType','frozen0',@(x)ismember(x, {'none', 'frozen', 'frozen0', 'random'}));
p.addParameter('responseStabilizationMilliseconds', 80, @isnumeric);
p.addParameter('responseExtinctionMilliseconds', 200, @isnumeric);
p.addParameter('computeResponses',true,@islogical);
p.addParameter('computeMosaic',false,@islogical);
p.addParameter('ramPercentageEmployed', 1.0, @isnumeric);
p.addParameter('parforWorkersNum', 12, @isnumeric);
% OPTICS OPTIONS
p.addParameter('pupilDiamMm',2,@isnumeric);
% MOSAIC OPTIONS
p.addParameter('freezeNoise',true, @islogical);
p.addParameter('osNoise', 'random', @(x)ismember(x, {'random', 'frozen', 'none'}));
p.addParameter('coneMosaicPacking', 'hex', @(x)ismember(x, {'hex', 'hexReg', 'rect'}));
p.addParameter('coneMosaicFOVDegs', 1.0, @isnumeric);
p.addParameter('integrationTime', 6.0/1000, @isnumeric);
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
p.addParameter('spatialPoolingKernelParams', struct(), @isstruct);
p.addParameter('performanceClassifier', 'mlpt', @(x)ismember(x, {'svm', 'svmSpaceTimeSeparable', 'svmGaussianRF', 'mlpt', 'mlpe', 'mlgtGaussianRF'}));
p.addParameter('performanceSignal', 'isomerizations', @(x)ismember(x, {'isomerizations', 'photocurrents'}));
p.addParameter('performanceClassifierTrainingSamples', [], @isnumeric);
p.addParameter('performanceEvidenceIntegrationTime', [], @isnumeric);
p.addParameter('summaryCurveMarkerType', 'o', @ischar);
p.addParameter('summaryCurveLegend', '', @ischar);
p.addParameter('plotSummaryCurve', false, @islogical);
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
stimulusDurationInSeconds = p.Results.stimulusDurationSecs;
if (strcmp(p.Results.stimulusTemporalEnvelope, 'Gaussian'))
    frameRate = 60;
    windowTauInSeconds = stimulusDurationInSeconds/3.0;
    stimulusSamplingIntervalInSeconds = 1/frameRate;
elseif (strcmp(p.Results.stimulusTemporalEnvelope, 'square'))
    frameRate = 200;
    windowTauInSeconds = nan; % square-wave
    stimulusSamplingIntervalInSeconds = 1/frameRate;
end

% Allow around 80 milliseconds for response to stabilize
responseStabilizationSeconds = ceil(p.Results.responseStabilizationMilliseconds/1000/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
% Allow around 200 milliseconds for response to return to 0
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

% Modify optics parameters
rParams.oiParams = modifyStructParams(rParams.oiParams, ...
            'pupilDiamMm', p.Results.pupilDiamMm);
        
% Modify mosaic parameters
rParams.mosaicParams = modifyStructParams(rParams.mosaicParams, ...
    'conePacking', p.Results.coneMosaicPacking, ...                       
    'fieldOfViewDegs', p.Results.coneMosaicFOVDegs, ... 
    'realisticSconeSubmosaic', true, ...             % if true, there will be no S-cones in the central 0.3 degs, and the S-cone lattice will be semiregular
    'integrationTimeInSeconds', p.Results.integrationTime, ...
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

% Adjust pca components num
pcaComponentsNum = p.Results.pcaComponentsNum;
if (strcmp(p.Results.performanceClassifier, 'svmSpaceTimeSeparable'))
    pcaComponentsNum = 4;
end

thresholdParams = modifyStructParams(thresholdParams, ...
    'method', p.Results.performanceClassifier, ...
    'STANDARDIZE', false, ...
    'standardizeSVMpredictors', false, ...
    'PCAComponents', pcaComponentsNum, ...
    'evidenceIntegrationTime', p.Results.performanceEvidenceIntegrationTime, ...
    'signalSource', p.Results.performanceSignal);

if (strcmp(thresholdParams.method, 'svm')) || ...
   (strcmp(thresholdParams.method, 'svmV1FilterBank')) || ...
   (strcmp(thresholdParams.method, 'svmGaussianRF'))
     thresholdParams = modifyStructParams(thresholdParams, ...
         'spatialPoolingKernelParams', p.Results.spatialPoolingKernelParams, ...
         'useRBFKernel', false);   % set to true to use a nonlinear (radial basis function) - based SVM
end

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
    baseLum = 12.5;
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
    if ((p.Results.visualizeResponses) || (p.Results.visualizeOuterSegmentFilters)) || (p.Results.visualizeSpatialScheme)
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
    if (p.Results.findPerformance) || (p.Results.visualizePerformance) || (p.Results.fitPsychometric) || (p.Results.visualizeSpatialScheme) || (p.Results.plotSummaryCurve)
        rParams.plotParams = modifyStructParams(rParams.plotParams, ...
            'axisFontSize', 12, ...
            'labelFontSize', 14, ...
            'lineWidth', 1.5);

        if (p.Results.findPerformance) || (p.Results.visualizePerformance)
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
        if (p.Results.fitPsychometric) && ((p.Results.findPerformance) || (p.Results.visualizePerformance)) || (p.Results.plotSummaryCurve)
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
    pupilAreaMM2 = pi * (p.Results.pupilDiamMm/2)^2;
    pedestalIlluminanceList = lumToTrolands(pedestalLuminanceList, pupilAreaMM2);
    hFig = figure(999); clf;
    set(hFig, 'Position', [10 10 450 450], 'Color', [1 1 1]);
    plot(log10(pedestalIlluminanceList), log10(irAmplitudes(:,1)), 'rs-', 'LineWidth', 1.5, 'MarkerSize', 14, 'MarkerFaceColor', [0.6 0.6 0.6]);
    hold on
    plot(log10(pedestalIlluminanceList), log10(irAmplitudes(:,2)), 'gs-', 'LineWidth', 1.5, 'MarkerSize', 14, 'MarkerFaceColor', [0.6 0.6 0.6]);
    plot(log10(pedestalIlluminanceList), log10(irAmplitudes(:,3)), 'bs-', 'LineWidth', 1.5, 'MarkerSize', 14, 'MarkerFaceColor', [0.6 0.6 0.6]);
    set(gca, 'FontSize', 14);
    set(gca, 'XLim', [log10(pedestalIlluminanceList(1)*0.95) log10(pedestalIlluminanceList(end)*1.05)]);
    grid on;
    box off;
    axis 'square'
    xlabel('log pedestal illuminance (Td)', 'FontWeight', 'bold');
    ylabel('log outer segment impulse response amplitude',  'FontWeight', 'bold');
    
    paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams};
    theProgram = mfilename;
    rwObject = IBIOColorDetectReadWriteBasic;
    data = 0;
    rwObject.write('os_ImpulseResponseSummary', data, paramsList, theProgram, ...
        'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
end


plotType = 'loglog';
%plotType = 'linear';

if (p.Results.plotSummaryCurve)
    detectionThresholdContrast = cell2mat(detectionThresholdContrast);
    thresholdLuminance = pedestalLuminanceList .* (detectionThresholdContrast);
    
    pupilAreaMM2 = pi * (p.Results.pupilDiamMm/2)^2;
    pedestalIlluminanceList = LumToTrolands(pedestalLuminanceList, pupilAreaMM2);
    thresholdIlluminance = LumToTrolands(thresholdLuminance, pupilAreaMM2);
    
    YLimContrast = [0.001 10];
    YLimIlluminance = [0.1 1000];
    XLimIlluminance = [1 10000];
    
    if (strcmp(plotType, 'loglog'))
        pedestalIlluminanceList = log10(pedestalIlluminanceList);
        detectionThresholdContrast = log10(detectionThresholdContrast);
        thresholdIlluminance = log10(thresholdIlluminance);
        XLimIlluminance = log10(XLimIlluminance);
        YLimContrast = log10(YLimContrast);
        contrastTicks = log10([0.01:0.01:0.2]);
        YLimIlluminance = log10(YLimIlluminance);
        IlluminanceTicks = log10(1:10);
    else 
        contrastTicks = [0.01:0.01:0.2];
        illuminanceTicks = 1:10;
    end
    
    
    GeislerData.x  = [-0.245    0.2914    0.7881    1.444     1.94      2.477     3.013   3.43    3.927    4.583];
    GeislerData.y1 = [-0.3487  -0.25    -0.072      0.2362    0.5395    0.9934    1.368   1.763   2.197    2.928];
    GeislerData.y2 = [0.2237    0.3421   0.4408     0.7961    1.013     1.428     1.783   2.237   2.711    3.283];
    GeislerData.y3 = [0.937     1.092    1.23       1.467     1.586     1.98      2.336   2.73    3.164    3.836];
    
    hFig = figure(1000); clf;
    set(hFig, 'Position', [10 10 450 900], 'Color', [1 1 1]);
    
    
    markerFaceColor1 = [0.7 0.7 0.7];
    markerEdgeColor1 = [0 0 0];
    markerFaceColor2 = [0.7 0.7 0.7];
    markerEdgeColor2 = [0 0 0];
    markerFaceColor3 = [0.7 0.7 0.7];
    markerEdgeColor3 = [0 0 0];
        
    subplot(2,1,1);
    hold on;
    if (strcmp(p.Results.summaryCurveMarkerType, 'o'))
        markerFaceColor2 = [0.9 0.9 0.9];
        markerEdgeColor2 = [0.5 0.5 0.5];
        markerFaceColor3 = [0.9 0.9 0.9];
        markerEdgeColor3 = [0.5 0.5 0.5]; 
    elseif (strcmp(p.Results.summaryCurveMarkerType, 's'))
        markerFaceColor1 = [0.9 0.9 0.9];
        markerEdgeColor1 = [0.5 0.5 0.5];
        markerFaceColor3 = [0.9 0.9 0.9];
        markerEdgeColor3 = [0.5 0.5 0.5];
    elseif (strcmp(p.Results.summaryCurveMarkerType, '^'))
        markerFaceColor1 = [0.9 0.9 0.9];
        markerEdgeColor1 = [0.5 0.5 0.5];
        markerFaceColor2 = [0.9 0.9 0.9];
        markerEdgeColor2 = [0.5 0.5 0.5];
    end
    if (strcmp(p.Results.summaryCurveMarkerType, 'o'))
        plot(GeislerData.x, GeislerData.y1, 'ko-', 'MarkerSize', 12, 'MarkerEdgeColor', markerEdgeColor1, 'MarkerFaceColor', markerFaceColor1, 'LineWidth', 1.5);
    end
    if (strcmp(p.Results.summaryCurveMarkerType, 's'))
        plot(GeislerData.x, GeislerData.y2, 'ks-', 'MarkerSize', 12, 'MarkerEdgeColor', markerEdgeColor2, 'MarkerFaceColor', markerFaceColor2, 'LineWidth', 1.5);
    end
    if (strcmp(p.Results.summaryCurveMarkerType, '^'))
        plot(GeislerData.x, GeislerData.y3, 'k^-', 'MarkerSize', 12, 'MarkerEdgeColor', markerEdgeColor3, 'MarkerFaceColor', markerFaceColor3, 'LineWidth', 1.5);
    end
    plot(pedestalIlluminanceList, thresholdIlluminance, 'r-', 'Marker', p.Results.summaryCurveMarkerType, 'LineWidth', 1.5, 'MarkerSize', 14, 'MarkerFaceColor', [1 0.8 0.8]);
    hold off
    axis 'equal'
    hL = legend({'Geisler (10'')', p.Results.summaryCurveLegend});
    set(hL, 'FontSize', 12, 'Location', 'NorthWest');
    set(gca, 'XLim', XLimIlluminance, 'YLim', YLimIlluminance, 'XTick', -5:5, 'YTick', -5:5);
    title(sprintf('%s (%s)',thresholdParams.signalSource, thresholdParams.method))
    if (strcmp(plotType, 'loglog'))
        xlabel('log adapting illuminance (Td)', 'FontWeight', 'bold');
        ylabel('log threshold test illuminance (Td)', 'FontWeight', 'bold');
    else
        xlabel('adapting illuminance (Td)', 'FontWeight', 'bold');
        ylabel('threshold test illuminance (Td)', 'FontWeight', 'bold');
    end
    set(gca, 'FontSize', 14);
    grid on;
    box on;
    
    
    subplot(2,1,2)
    plot(pedestalIlluminanceList, detectionThresholdContrast, 'r-', 'Marker', p.Results.summaryCurveMarkerType, 'LineWidth', 1.5, 'MarkerSize', 14, 'MarkerFaceColor', [1 0.8 0.8]);
    axis 'equal'
    set(gca, 'XLim', XLimIlluminance, 'YLim', YLimContrast, 'XTick', -5:5, 'YTick', -5:5);
    if (strcmp(plotType, 'loglog'))
        xlabel('log adapting illuminance (Td)', 'FontWeight', 'bold');
        ylabel('log threshold test contrast', 'FontWeight', 'bold');
    else
        xlabel('adapting illuminance (Td)', 'FontWeight', 'bold');
        ylabel('threshold test contrast', 'FontWeight', 'bold');
    end
    title(sprintf('%s (%s)',thresholdParams.signalSource, thresholdParams.method));
    set(gca, 'FontSize', 14);
    grid on;
    box on;
    
    
    
    
    drawnow

    
    paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams, thresholdParams};
    theProgram = mfilename;
    rwObject = IBIOColorDetectReadWriteBasic;
    data = 0;
    rwObject.write('c_TCIcurveSummary', data, paramsList, theProgram, ...
        'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
end
       
end

