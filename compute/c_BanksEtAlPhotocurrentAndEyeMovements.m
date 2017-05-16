function varargout = c_BanksEtAlPhotocurrentAndEyeMovements(varargin)

    
%% Parse input
p = inputParser;
% ----- cBanksEtAl params -----
p.addParameter('employStandardHostComputerResources', false, @islogical);
p.addParameter('useScratchTopLevelDirName', false, @islogical);
p.addParameter('nTrainingSamples',500,@isnumeric);
p.addParameter('cyclesPerDegree',[3 5 10 20 40 50],@isnumeric); 
p.addParameter('spatialPhaseDegs',0,@isnumeric);
p.addParameter('luminances',[3.4 34 340],@isnumeric);
p.addParameter('pupilDiamMm',2,@isnumeric);
p.addParameter('blur',true,@islogical);
p.addParameter('opticsModel','WvfHuman',@ischar);
p.addParameter('innerSegmentSizeMicrons',sizeForSquareApertureFromDiameterForCircularAperture(3.0), @isnumeric);   % 3 microns = 0.6 min arc for 300 microns/deg in human retina
p.addParameter('apertureBlur', false, @islogical);
p.addParameter('coneSpacingMicrons', 3.0, @isnumeric);
p.addParameter('mosaicRotationDegs', 0, @isnumeric);
p.addParameter('coneDarkNoiseRate',[0 0 0], @isnumeric);
p.addParameter('LMSRatio',[0.67 0.33 0],@isnumeric);
p.addParameter('conePacking', 'hexReg',@ischar);                 
p.addParameter('imagePixels',400,@isnumeric);
p.addParameter('wavelengths',[380 4 780],@isnumeric);
p.addParameter('nContrastsPerDirection',20,@isnumeric);
p.addParameter('lowContrast',0.0001,@isnumeric);
p.addParameter('highContrast',0.1,@isnumeric);
p.addParameter('contrastScale','log',@ischar);
p.addParameter('computeResponses',true,@islogical);
p.addParameter('computePhotocurrentResponseInstances', true, @islogical);
p.addParameter('findPerformance',true,@islogical);
p.addParameter('thresholdMethod','mlpt',@ischar);
p.addParameter('thresholdPCA',60,@isnumeric);
p.addParameter('fitPsychometric',true,@islogical);
p.addParameter('thresholdCriterionFraction',0.75,@isnumeric);
p.addParameter('generatePlots',true,@islogical);
p.addParameter('visualizeResponses',false,@islogical);
p.addParameter('visualizedConditionIndices', [], @isnumeric);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.addParameter('plotPsychometric',true,@islogical);
p.addParameter('plotCSF',true,@islogical);
p.addParameter('freezeNoise',true,@islogical);

% --- additional params -----
% RESPONSE COMPUTATION OPTIONS
p.addParameter('emPathType','frozen0',@(x)ismember(x, {'none', 'frozen', 'frozen0', 'random'}));
p.addParameter('centeredEMPaths',false, @islogical); 
p.addParameter('responseStabilizationMilliseconds', 80, @isnumeric);
p.addParameter('responseExtinctionMilliseconds', 200, @isnumeric);
p.addParameter('computeMosaic',false,@islogical);
p.addParameter('ramPercentageEmployed', 1.0, @isnumeric);
p.addParameter('parforWorkersNum', 20, @isnumeric);
p.addParameter('parforWorkersNumForClassification', 6, @isnumeric);

% MOSAIC OPTIONS
p.addParameter('integrationTime', 5.0/1000, @isnumeric);

% DIAGNOSTIC OPTIONS
p.addParameter('displayTrialBlockPartitionDiagnostics', true, @islogical);
p.addParameter('displayResponseComputationProgress', false, @islogical);

% RESPONSE MAP VISUALIZATION OPTIONS
p.addParameter('visualizeMosaic',true, @islogical); 
p.addParameter('visualizeSpatialScheme', false, @islogical);
p.addParameter('visualizeTransformedSignals',false, @islogical);
p.addParameter('visualizationFormat', 'montage', @ischar);
p.addParameter('visualizePerformance', false, @islogical);

% PERFORMANCE COMPUTATION OPTIONS
p.addParameter('spatialPoolingKernelParams', struct(), @isstruct);
p.addParameter('useRBFSVMKernel', false, @islogical);
p.addParameter('performanceClassifier', 'mlpt', @(x)ismember(x, {'svm', 'svmSpaceTimeSeparable', 'svmGaussianRF', 'svmV1FilterBank', 'mlpt', 'mlpe', 'mlgtGaussianRF'}));
p.addParameter('performanceSignal', 'isomerizations', @(x)ismember(x, {'isomerizations', 'photocurrents'}));
p.addParameter('performanceTrialsUsed', [], @isnumeric);

% DELETE RESPONSE INSTANCES
p.addParameter('deleteResponseInstances', false, @islogical);

p.parse(varargin{:});



%% Set the default output
varargout = {};
varargout{1} = [];
varargout{2} = [];
varargout{3} = [];
   
%% Make sure we really want to delete the response files
if (p.Results.deleteResponseInstances)
    proceed = input('Responses will be deleted. Are you sure ? [y = YES]', 's');
    if (~strcmp(proceed, 'y'))
        return;
    end
end

%% Get the parameters we need
%
% Start with default
rParams = responseParamsGenerate;

%% Set the  topLevelDir name
if (~p.Results.useScratchTopLevelDirName)
    rParams.topLevelDirParams.name = mfilename;
end

%% Loop over spatial frequency
for ll = 1:length(p.Results.luminances)
    for cc = 1:length(p.Results.cyclesPerDegree)
        
        % Modify spatial params
        cyclesPerDegree = p.Results.cyclesPerDegree(cc);
        gaussianFWHMDegs = 3.75*(1/cyclesPerDegree);
        fieldOfViewDegs = 2.1*gaussianFWHMDegs;
        rParams.spatialParams = modifyStructParams(rParams.spatialParams, ...
            'windowType', 'halfcos', ...
            'cyclesPerDegree', cyclesPerDegree, ...
            'ph', (pi/180)*p.Results.spatialPhaseDegs, ...
            'gaussianFWHMDegs', gaussianFWHMDegs, ...
            'fieldOfViewDegs', fieldOfViewDegs, ...
            'row', p.Results.imagePixels, ...
            'col', p.Results.imagePixels);
        
        % Set wavelength sampling
        rParams.colorModulationParams.startWl = p.Results.wavelengths(1);
        rParams.colorModulationParams.deltaWl = p.Results.wavelengths(2);
        rParams.colorModulationParams.endWl = p.Results.wavelengths(3);
        
        % Blur
        %
        % Banks et al. used a 2mm artificial pupil
        rParams.oiParams = modifyStructParams(rParams.oiParams, ...
        	'blur', p.Results.blur, ...
            'pupilDiamMm', p.Results.pupilDiamMm, ...
            'opticsModel', p.Results.opticsModel);
        
        % Modify background params
        baseLum = 50;
        theLum = p.Results.luminances(ll);
        rParams.backgroundParams = modifyStructParams(rParams.backgroundParams, ...
        	'backgroundxyY', [0.33 0.33 baseLum]',...
        	'monitorFile', 'CRT-MODEL', ...
        	'leakageLum', 1.0, ...
        	'lumFactor', theLum/baseLum);
        
        % Temporal params
        % Their stimulus intervals were 100 msec each.
        stimulusDurationInSeconds = 100/1000;
        
        % In the absence of info from the Banks et al paper, assume 50 Hz refresh rate
        frameRate = 50;
        windowTauInSeconds = nan; % square-wave
        stimulusSamplingIntervalInSeconds = 1/frameRate;
    
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

        % Modify mosaic parameters
        rParams.mosaicParams = modifyStructParams(rParams.mosaicParams, ...
            'conePacking', p.Results.conePacking, ...                       
            'fieldOfViewDegs', rParams.spatialParams.fieldOfViewDegs, ... 
            'realisticSconeSubmosaic', true, ...             % if true, there will be no S-cones in the central 0.3 degs, and the S-cone lattice will be semiregular
            'innerSegmentSizeMicrons',p.Results.innerSegmentSizeMicrons, ...
            'coneSpacingMicrons', p.Results.coneSpacingMicrons, ...
            'apertureBlur',p.Results.apertureBlur, ...
            'mosaicRotationDegs',p.Results.mosaicRotationDegs,...
            'coneDarkNoiseRate',p.Results.coneDarkNoiseRate,...
            'LMSRatio',p.Results.LMSRatio,...
            'integrationTimeInSeconds', p.Results.integrationTime, ...
            'isomerizationNoise', 'random',...               % select from {'random', 'frozen', 'none'}
            'osNoise', 'random', ...                % select from {'random', 'frozen', 'none'}
            'osModel', 'Linear');


        % Make sure mosaic noise parameters match the frozen noise flag passed in
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
        testDirectionParams = instanceParamsGenerate;
        testDirectionParams = modifyStructParams(testDirectionParams, ...
            'trialsNum', p.Results.nTrainingSamples, ...
        	'startAngle', 45, ...
        	'deltaAngle', 90, ...
        	'nAngles', 1, ...
            'nContrastsPerDirection', p.Results.nContrastsPerDirection, ...
            'lowContrast', p.Results.lowContrast, ...
        	'highContrast', p.Results.highContrast, ...
        	'contrastScale', p.Results.contrastScale ...                       % choose between 'linear' and 'log'
            );
        
        % Parameters related to how we find thresholds from responses
        thresholdParams = thresholdParamsGenerate;
        thresholdParams = modifyStructParams(thresholdParams, ...
        	'criterionFraction', p.Results.thresholdCriterionFraction, ...
        	'method', p.Results.performanceClassifier, ...
            'STANDARDIZE', false, ...                   % Standardize data for PCA
            'standardizeSVMpredictors', false, ...       % Standardize data for SVM
            'useRBFKernel', p.Results.useRBFSVMKernel, ...
        	'PCAComponents', p.Results.thresholdPCA, ...
            'signalSource', p.Results.performanceSignal ...
            );
        if (strcmp(thresholdParams.method, 'svm')) || ...
            (strcmp(thresholdParams.method, 'svmV1FilterBank')) || ...
            (strcmp(thresholdParams.method, 'svmGaussianRF'))
            thresholdParams = modifyStructParams(thresholdParams, ...
                'spatialPoolingKernelParams', p.Results.spatialPoolingKernelParams ...
                );  
        end

        if (isempty(p.Results.performanceTrialsUsed))
            thresholdParams.trialsUsed = p.Results.nTrainingSamples;
        else
            thresholdParams.trialsUsed = p.Results.performanceTrialsUsed;
        end

        %% Compute response instances
        if (p.Results.computeResponses) || (p.Results.visualizeMosaic) 
           t_coneCurrentEyeMovementsResponseInstances(...
               'rParams',rParams,...
               'testDirectionParams',testDirectionParams,...
               'centeredEMPaths',p.Results.centeredEMPaths, ...
               'compute',p.Results.computeResponses, ...
               'computePhotocurrentResponseInstances', p.Results.computePhotocurrentResponseInstances, ...
               'computeMosaic', (ll == 1) && (p.Results.computeMosaic), ...   % only compute the mosaic for the first luminance level 
               'visualizeMosaic', p.Results.visualizeMosaic, ...
               'visualizeResponses',false, ...
               'visualizeSpatialScheme', p.Results.visualizeSpatialScheme, ...
               'generatePlots',p.Results.generatePlots,...
               'freezeNoise',p.Results.freezeNoise,...
               'ramPercentageEmployed', p.Results.ramPercentageEmployed, ...
               'parforWorkersNum', p.Results.parforWorkersNum, ...  % no more than these many workers
               'displayTrialBlockPartitionDiagnostics', p.Results.displayTrialBlockPartitionDiagnostics, ...
               'employStandardHostComputerResources', p.Results.employStandardHostComputerResources, ...
               'workerID', find(p.Results.displayResponseComputationProgress) ...
            );
        end

        %% Visualize response instances
        if (p.Results.visualizeResponses)
            [~, ~, ~, noiseFreeResponses{cc,ll}] = t_coneCurrentEyeMovementsResponseInstances(...
              'rParams',rParams,...
              'testDirectionParams',testDirectionParams,...
              'freezeNoise', p.Results.freezeNoise, ...
              'compute', false, ...
              'computeMosaic', false, ... 
              'visualizeMosaic', p.Results.visualizeMosaic, ...
              'visualizedResponseNormalization', p.Results.visualizedResponseNormalization, ...
              'visualizationFormat', p.Results.visualizationFormat, ...
              'generatePlots', p.Results.generatePlots, ...
              'visualizeResponses', p.Results.visualizeResponses, ...
              'visualizedConditionIndices', p.Results.visualizedConditionIndices, ...
              'visualizeOuterSegmentFilters', false ...
            );
        end
            
        %% Find performance
        if (p.Results.findPerformance) || (p.Results.visualizePerformance)
            t_colorDetectFindPerformance(...
                'rParams',rParams, ...
                'testDirectionParams',testDirectionParams, ...
                'thresholdParams',thresholdParams, ...
                'compute',p.Results.findPerformance, ...
                'parforWorkersNum', p.Results.parforWorkersNumForClassification, ...
                'visualizeSpatialScheme', p.Results.visualizeSpatialScheme, ...
                'visualizeTransformedSignals', p.Results.visualizeTransformedSignals, ...
                'plotSvmBoundary',false,...
                'plotPsychometric',false,...
                'freezeNoise',p.Results.freezeNoise ...
            );
        
        
            %% Fit psychometric functions
            if (p.Results.fitPsychometric)
                banksEtAlReplicate.cyclesPerDegree(ll,cc) = p.Results.cyclesPerDegree(cc);
                banksEtAlReplicate.mlptThresholds(ll,cc) = ...
                    t_plotDetectThresholdsOnLMPlane(...
                        'rParams',rParams, ...
                        'instanceParams',testDirectionParams, ...
                        'thresholdParams',thresholdParams, ...
                        'plotPsychometric', p.Results.generatePlots & p.Results.plotPsychometric, ...
                        'plotEllipse',false);
                %close all;
            end
        end
        
        %% Delete response instances
        if (p.Results.deleteResponseInstances)
            t_coneCurrentEyeMovementsResponseInstances(...
               'rParams',rParams,...
               'testDirectionParams',testDirectionParams,...
               'centeredEMPaths',p.Results.centeredEMPaths, ...
               'freezeNoise',p.Results.freezeNoise,...
               'compute',false, ...
               'computePhotocurrentResponseInstances', false, ...
               'computeMosaic', false, ...   
               'visualizeMosaic', false, ...
               'visualizeResponses',false, ...
               'visualizeSpatialScheme', false, ...
               'generatePlots', false,...
                'delete', true);
        end
        
        rParamsAllConds{cc,ll} = rParams;
    end % cc
end % ll

if (p.Results.visualizeResponses)
    varargout{3} = noiseFreeResponses;
end

% Write out the data
if (p.Results.fitPsychometric) && ((p.Results.findPerformance) || (p.Results.visualizePerformance))
    
    fprintf('Writing performance data ... ');
    paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams, thresholdParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    rwObject.write('banksEtAlReplicatePhotocurrentAndEyeMovements',banksEtAlReplicate,paramsList,writeProgram);
    fprintf('done\n');


    %% Get performance data
    fprintf('Reading performance data ...');
    paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams, thresholdParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    banksEtAlReplicate = rwObject.read('banksEtAlReplicatePhotocurrentAndEyeMovements',paramsList,writeProgram);
    fprintf('done\n');

    % Return the performance data
    varargout{1} = banksEtAlReplicate;
    varargout{2} = rParamsAllConds;
    
    %% Make a plot of CSFs
    %
    % The way the plot is coded counts on the test contrasts never changing
    % across the conditions, which we could explicitly check for here.
    if (p.Results.generatePlots && p.Results.plotCSF)  
        hFig = figure; clf; hold on
        fontBump = 4;
        markerBump = -4;
        lineBump = -1;
        set(gcf,'Position',[100 100 450 650]);
        set(gca,'FontSize', rParams.plotParams.axisFontSize+fontBump);
        theColors = ['r' 'g' 'b'];
        legendStr = cell(length(p.Results.luminances),1);
        for ll = 1:length(p.Results.luminances)
            theColorIndex = 3; % rem(ll,length(theColors)) + 1;
            plot(banksEtAlReplicate.cyclesPerDegree(ll,:),1./([banksEtAlReplicate.mlptThresholds(ll,:).thresholdContrasts]*banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1)), ...
                [theColors(theColorIndex) 'o'],'MarkerSize',rParams.plotParams.markerSize+markerBump,'MarkerFaceColor',theColors(theColorIndex),'LineWidth',rParams.plotParams.lineWidth);  
            legendStr{ll} = sprintf('%0.1f cd/m2',p.Results.luminances(ll));
        end

        % Add Banks et al. data to the figure
        % Slide by eye to match our current calculations, and add as
        % appropriate so as not to clutter the figure.  A bit klugy as we add
        % more and more conditions, I fear.
        % 
        % Also add unshifted version for reference
        banksFactor = 1;
        [A,B,C,D,E] = LoadDigitizedBanksFigure2;
        if (~rParams.oiParams.blur && ~rParams.mosaicParams.apertureBlur)
            plot(A(:,1),A(:,2),'k:','LineWidth',0.5);
            plot(A(:,1),A(:,2)*banksFactor,'r-','LineWidth',rParams.plotParams.lineWidth+lineBump);
        elseif (~rParams.oiParams.blur)
            plot(B(:,1),B(:,2),'k:','LineWidth',0.5);
            plot(B(:,1),B(:,2)*banksFactor,'r','LineWidth',rParams.plotParams.lineWidth+lineBump);
        else
            plot(C(:,1),C(:,2),'k:','LineWidth',0.5);
            plot(C(:,1),C(:,2)*banksFactor,'r-','LineWidth',rParams.plotParams.lineWidth+lineBump);
            plot(D(:,1),D(:,2)*banksFactor,'b-','LineWidth',rParams.plotParams.lineWidth+lineBump);
            plot(E(:,1),E(:,2)*banksFactor,'g-','LineWidth',rParams.plotParams.lineWidth+lineBump);
        end

        set(gca,'XScale','log');
        set(gca,'YScale','log');
        xlabel('Log10 Spatial Frequency (cpd)', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
        ylabel('Log10 Contrast Sensitivity', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
        xlim([1 100]); ylim([1 10000]);
        legend(legendStr,'Location','NorthEast','FontSize',rParams.plotParams.labelFontSize+fontBump);
        box off; grid on
        titleStr1 = 'Computational Observer CSF';
        titleStr2 = sprintf('Blur: %d, Aperture Blur: %d',p.Results.blur, p.Results.apertureBlur);
        title({titleStr1 ; titleStr2});

        % Write out the figure
        rwObject.write('banksEtAlReplicatePhotocurrentAndEyeMovements',hFig,paramsList,writeProgram,'Type','figure')
    end
end


end

