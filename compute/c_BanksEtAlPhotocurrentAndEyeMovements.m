function varargout = c_BanksEtAlPhotocurrentAndEyeMovements(varargin)

    
%% Parse input
p = inputParser;

% ----- cBanksEtAl params -----
p.addParameter('employStandardHostComputerResources', false, @islogical);
p.addParameter('useScratchTopLevelDirName', false, @islogical);
p.addParameter('nTrainingSamples',500,@isnumeric);
p.addParameter('cyclesPerDegree',[2.5 5 10 20 40 60],@isnumeric); 
p.addParameter('spatialPhaseDegs',0,@isnumeric);
p.addParameter('luminances',[3.4 34 340],@isnumeric);

% Optics params
p.addParameter('pupilDiamMm',2,@isnumeric);
p.addParameter('blur',true,@islogical);
p.addParameter('opticsModel','WvfHuman',@ischar);
p.addParameter('wavefrontSpatialSamples', 201, @isnumeric);      %% * * * NEW * * * 
p.addParameter('minimumOpticalImagefieldOfViewDegs', 1.0, @isnumeric);  %% * * * NEW * * * 

% Mosaic params
p.addParameter('integrationTime', 5.0/1000, @isnumeric);
p.addParameter('innerSegmentSizeMicrons',3.0, @isnumeric);   % 3 microns = 0.6 min arc for 300 microns/deg in human retina
p.addParameter('apertureBlur',true, @islogical);
p.addParameter('coneSpacingMicrons', 3.0, @isnumeric);
p.addParameter('mosaicRotationDegs', 0, @isnumeric);
p.addParameter('coneDarkNoiseRate',[0 0 0], @isnumeric);
p.addParameter('LMSRatio',[0.67 0.33 0],@isnumeric);
p.addParameter('conePacking', 'hexReg',@ischar);
p.addParameter('eccBasedConeQuantalEfficiency', false, @islogical);
p.addParameter('freezeNoise',true,@islogical);

% HEX MOSAIC OPTIONS
p.addParameter('sConeMinDistanceFactor', 3.0, @isnumeric); % min distance between neighboring S-cones = f * local cone separation - to make the S-cone lattice semi-regular
p.addParameter('sConeFreeRadiusMicrons', 45, @isnumeric);
p.addParameter('latticeAdjustmentPositionalToleranceF', 0.01, @isnumeric);
p.addParameter('latticeAdjustmentDelaunayToleranceF', 0.001, @isnumeric);
p.addParameter('maxGridAdjustmentIterations', Inf, @isnumeric);  % * * * NEW * * * 
p.addParameter('marginF', [], @isnumeric);     
p.addParameter('resamplingFactor', 9, @isnumeric);              

% Stimulus params
p.addParameter('imagePixels',400,@isnumeric);
p.addParameter('wavelengths',[380 5 780],@isnumeric);
p.addParameter('coneContrastDirection', 'L+M', @ischar);
p.addParameter('stimulusDurationInSeconds', 100/1000, @isnumeric);   
p.addParameter('nContrastsPerDirection',20,@isnumeric);
p.addParameter('lowContrast',0.0001,@isnumeric);
p.addParameter('highContrast',0.1,@isnumeric);
p.addParameter('contrastScale','log',@ischar);
p.addParameter('frameRate',10,@isnumeric);

% Response dynamics
p.addParameter('emPathType','frozen0',@(x)ismember(x, {'none', 'frozen', 'frozen0', 'random', 'randomNoSaccades'}));
p.addParameter('centeredEMPaths',false, @islogical); 
p.addParameter('responseStabilizationMilliseconds', 80, @isnumeric);
p.addParameter('responseExtinctionMilliseconds', 200, @isnumeric);

% PERFORMANCE COMPUTATION OPTIONS
p.addParameter('spatialPoolingKernelParams', struct(), @isstruct);
p.addParameter('useRBFSVMKernel', false, @islogical);
p.addParameter('performanceClassifier', 'mlpt', @(x)ismember(x, {'svm', 'svmSpaceTimeSeparable', 'svmGaussianRF', 'svmV1FilterBank', 'svmV1FilterEnsemble', 'mlpt', 'mlpe', 'mlgtGaussianRF'}));
p.addParameter('performanceSignal', 'isomerizations', @(x)ismember(x, {'isomerizations', 'photocurrents'}));
p.addParameter('performanceTrialsUsed', [], @isnumeric);
p.addParameter('thresholdCriterionFraction',0.7071,@isnumeric);

% RESPONSE COMPUTATION OPTIONS
p.addParameter('ramPercentageEmployed', 1.0, @isnumeric); 
p.addParameter('parforWorkersNum', 20, @isnumeric);
p.addParameter('parforWorkersNumForClassification', 6, @isnumeric);

% What to compute
p.addParameter('computeOptics',true,@islogical);
p.addParameter('computeMosaic',false,@islogical);
p.addParameter('computeResponses',true,@islogical);
p.addParameter('computePhotocurrentResponseInstances', true, @islogical);
p.addParameter('findPerformance',true,@islogical);
p.addParameter('thresholdMethod','mlpt',@ischar);
p.addParameter('thresholdPCA',60,@isnumeric);
p.addParameter('fitPsychometric',true,@islogical);


% What to visualize
p.addParameter('generatePlots',true,@islogical);
p.addParameter('visualizeDisplay', false, @islogical);
p.addParameter('visualizeOIsequence', false, @islogical);
p.addParameter('visualizeStimulusAndOpticalImage', false, @islogical);
p.addParameter('visualizeOptics',false, @islogical);
p.addParameter('visualizeSpatialScheme', false, @islogical);
p.addParameter('visualizeSpatialPoolingScheme', false, @islogical);
p.addParameter('visualizeMosaic',false, @islogical); 
p.addParameter('visualizeMosaicWithFirstEMpath', false, @islogical);
p.addParameter('visualizeResponses',false,@islogical);
p.addParameter('visualizeKernelTransformedSignals',false, @islogical);
p.addParameter('visualizationFormat', 'montage', @ischar);
p.addParameter('visualizePerformance', false, @islogical);
p.addParameter('visualizedConditionIndices', [], @isnumeric);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.addParameter('plotPsychometric',true,@islogical);
p.addParameter('plotCSF',true,@islogical);

% DIAGNOSTIC OPTIONS
p.addParameter('displayTrialBlockPartitionDiagnostics', false, @islogical);
p.addParameter('displayResponseComputationProgress', false, @islogical);

% DELETE RESPONSE INSTANCES
p.addParameter('deleteResponseInstances', false, @islogical);

p.parse(varargin{:});

% Take a snapshot of the current IBIOColorDetect deployment
IBIOColorDetectConfig = tbUseProject('IBIOColorDetect');
IBIOColorDetectSnapshot = tbDeploymentSnapshot(IBIOColorDetectConfig);


% Following does not work, so we will write the stuct out ourselves
%tbWriteConfig(IBIOColorDetectShapshot, 'configPath', 'IBIOColorDetect_snapshot.json');

%% Set the default output
varargout = {};
varargout{1} = [];  % banksEtAlReplicate
varargout{2} = [];  % rParamsAllConds
varargout{3} = [];  % noise-free responses
varargout{4} = [];  % the optical image optics
varargout{5} = [];  % the cone mosaic
varargout{6} = {};  % the psychometric functions
varargout{7} = {};  % the figure data (CSF curves)

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

% Set wavelength sampling stimulus params
rParams = updateSpectralParams(rParams, p.Results);

% Update temporal stimulus params
rParams = updateTemporalParams(rParams, p.Results);

% Update eye movement params
rParams = updateEyeMovementParams(rParams, p.Results);

% Generate test direction params;
testDirectionParams = generateTestDirectionParams(p.Results);

% Generate threshold params
thresholdParams = generateThresholdParams(p.Results);

%% Loop over luminance and spatial frequency
for ll = 1:length(p.Results.luminances)
    for cc = 1:length(p.Results.cyclesPerDegree)
        
        % Modify spatial stimulus params
        rParams = updateSpatialParams(rParams, p.Results, p.Results.cyclesPerDegree(cc));
          
        % Modify background stimulus params
        rParams = updateBackgroundParams(rParams, p.Results, p.Results.luminances(ll));
        
        % Modify mosaic params
        rParams = updateMosaicParams(rParams, p.Results, rParams.spatialParams.fieldOfViewDegs);
        
        % Modify the optical image params
        rParams = updateOpticalImageParams(rParams, p.Results, rParams.spatialParams.fieldOfViewDegs);
    
        %% Compute response instances
        if (p.Results.computeResponses) || (p.Results.visualizeMosaic) || (p.Results.visualizeOptics) || (p.Results.computeMosaic)
           [validationData, extraData, ~,~, theOI, theMosaic, thePupilFunction] = t_coneCurrentEyeMovementsResponseInstances(...
               'rParams',rParams,...
               'testDirectionParams',testDirectionParams,...
               'centeredEMPaths',p.Results.centeredEMPaths, ...
               'compute',p.Results.computeResponses, ...
               'computePhotocurrentResponseInstances', p.Results.computePhotocurrentResponseInstances, ...
               'computeMosaic', (ll == 1) && (p.Results.computeMosaic), ...   % only compute the mosaic for the first luminance level
               'computeOptics', (ll == 1) && (p.Results.computeOptics), ...   % only compute the mosaic for the first luminance level 
               'visualizeMosaic', p.Results.visualizeMosaic, ...
               'visualizeMosaicWithFirstEMpath', p.Results.visualizeMosaicWithFirstEMpath, ...
               'visualizeDisplay', p.Results.visualizeDisplay, ...
               'visualizeOptics', p.Results.visualizeOptics, ...
               'visualizeStimulusAndOpticalImage', p.Results.visualizeStimulusAndOpticalImage, ...
               'visualizeOIsequence', p.Results.visualizeOIsequence, ...
               'visualizeResponses',false, ...
               'visualizeSpatialScheme', p.Results.visualizeSpatialScheme, ...
               'generatePlots',p.Results.generatePlots,...
               'freezeNoise',p.Results.freezeNoise,...
               'ramPercentageEmployed', p.Results.ramPercentageEmployed, ...
               'parforWorkersNum', p.Results.parforWorkersNum, ...  % no more than these many workers
               'displayTrialBlockPartitionDiagnostics', p.Results.displayTrialBlockPartitionDiagnostics, ...
               'employStandardHostComputerResources', p.Results.employStandardHostComputerResources, ...
               'workerID', find(p.Results.displayResponseComputationProgress), ...
               'IBIOColorDetectSnapshot', IBIOColorDetectSnapshot ...
            );
        
            %
            % Get the OI and the mosaic for the current spatial frequency
            varargout{4}.theOIs{cc} = theOI;
            varargout{5}.theMosaics{cc} = theMosaic;
            %end
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
              'visualizeOptics', p.Results.visualizeOptics, ...
              'visualizedResponseNormalization', p.Results.visualizedResponseNormalization, ...
              'visualizationFormat', p.Results.visualizationFormat, ...
              'generatePlots', p.Results.generatePlots, ...
              'visualizeResponses', p.Results.visualizeResponses, ...
              'visualizeMosaicWithFirstEMpath',false, ...
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
                'visualizeSpatialScheme', p.Results.visualizeSpatialPoolingScheme, ...
                'visualizeKernelTransformedSignals', p.Results.visualizeKernelTransformedSignals, ...
                'visualizeVarianceExplained', true, ...
                'plotSvmBoundary',false,...
                'plotPsychometric',false,...
                'freezeNoise',p.Results.freezeNoise, ...
                'IBIOColorDetectSnapshot', IBIOColorDetectSnapshot ...
            );
        
        
            %% Fit psychometric functions
            if (p.Results.fitPsychometric)
                banksEtAlReplicate.cyclesPerDegree(ll,cc) = p.Results.cyclesPerDegree(cc);
                [banksEtAlReplicate.mlptThresholds(ll,cc), psychometricFunctions{ll,cc}] = ...
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
               'visualizeOptics', false, ...
               'visualizeResponses',false, ...
               'visualizeSpatialScheme', false, ...
               'generatePlots', false,...
                'delete', true);
        end
        
        rParamsAllConds{cc,ll} = rParams;
    end % cc
end % ll

if ((p.Results.findPerformance) || (p.Results.visualizePerformance)) && (p.Results.fitPsychometric)
    varargout{6} = psychometricFunctions;
end

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
        figData.banksEtAlReplicate = banksEtAlReplicate;
        
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
            figData.A = A;
        elseif (~rParams.oiParams.blur)
            plot(B(:,1),B(:,2),'k:','LineWidth',0.5);
            plot(B(:,1),B(:,2)*banksFactor,'r','LineWidth',rParams.plotParams.lineWidth+lineBump);
            figData.B = B;
        else
            plot(C(:,1),C(:,2),'k:','LineWidth',0.5);
            plot(C(:,1),C(:,2)*banksFactor,'r-','LineWidth',rParams.plotParams.lineWidth+lineBump);
            plot(D(:,1),D(:,2)*banksFactor,'b-','LineWidth',rParams.plotParams.lineWidth+lineBump);
            plot(E(:,1),E(:,2)*banksFactor,'g-','LineWidth',rParams.plotParams.lineWidth+lineBump);
            figData.C = C;
            figData.D = D;
            figData.E = E;
        end

        set(gca,'XScale','log');
        set(gca,'YScale','log');
        xlabel('Log10 Spatial Frequency (cpd)', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
        ylabel('Log10 Contrast Sensitivity', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
        xlim([1 100]); ylim([1 12000]);
        legend(legendStr,'Location','NorthEast','FontSize',rParams.plotParams.labelFontSize+fontBump);
        box off; grid on
        titleStr1 = 'Computational Observer CSF';
        titleStr2 = sprintf('Blur: %d, Aperture Blur: %d',p.Results.blur, p.Results.apertureBlur);
        title({titleStr1 ; titleStr2});
        
        % Return the figure data
        varargout{7} = figData;
        % Write out the figure
        rwObject.write('banksEtAlReplicatePhotocurrentAndEyeMovements',hFig,paramsList,writeProgram,'Type','figure')
    end
end
end

% Method to update the spatial stimulus params
function rParams = updateSpatialParams(rParams, userParams, cyclesPerDegree)

    gaussianFWHMDegs = 3.75*(1/cyclesPerDegree);
    fieldOfViewDegs = 2.1*gaussianFWHMDegs;
    rParams.spatialParams = modifyStructParams(rParams.spatialParams, ...
            'windowType', 'halfcos', ...
            'cyclesPerDegree', cyclesPerDegree, ...
            'ph', (pi/180)*userParams.spatialPhaseDegs, ...
            'gaussianFWHMDegs', gaussianFWHMDegs, ...
            'fieldOfViewDegs', fieldOfViewDegs, ...
            'row', userParams.imagePixels, ...
            'col', userParams.imagePixels);
end

% Method to update the spectral stimulus params
function rParams = updateSpectralParams(rParams, userParams)
    rParams.colorModulationParams.startWl = userParams.wavelengths(1);
    rParams.colorModulationParams.deltaWl = userParams.wavelengths(2);
    rParams.colorModulationParams.endWl = userParams.wavelengths(3);
end

% Method to update the background stimulus params
function rParams = updateBackgroundParams(rParams, userParams, luminance)
    baseLum = 40;
    rParams.backgroundParams = modifyStructParams(rParams.backgroundParams, ...
        	'backgroundxyY', [0.33 0.33 baseLum]',...
        	'monitorFile', 'CRT-MODEL', ...
        	'leakageLum', 1.0, ...
        	'lumFactor', luminance/baseLum);
        
end

function rParams = updateTemporalParams(rParams, userParams)
    % In the absence of info from the Banks et al paper, assume 10 Hz
    % refresh rate (to accelerate computations)
    frameRate = userParams.frameRate;
    windowTauInSeconds = nan; % square-wave
    stimulusSamplingIntervalInSeconds = 1/frameRate;
    
    % Allow some milliseconds for response to stabilize
    responseStabilizationSeconds = ceil(userParams.responseStabilizationMilliseconds/1000/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
    % Allow some milliseconds for response to return to 0
    responseExtinctionSeconds = ceil(userParams.responseExtinctionMilliseconds/1000/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
    secondsToInclude = responseStabilizationSeconds+userParams.stimulusDurationInSeconds+responseExtinctionSeconds;   
        
    rParams.temporalParams = modifyStructParams(rParams.temporalParams, ...
            'frameRate', frameRate, ...
            'windowTauInSeconds', windowTauInSeconds, ...
            'stimulusSamplingIntervalInSeconds', stimulusSamplingIntervalInSeconds, ...
            'stimulusDurationInSeconds', userParams.stimulusDurationInSeconds, ...
            'secondsToInclude', secondsToInclude, ...
            'secondsForResponseStabilization', responseStabilizationSeconds, ...
            'secondsForResponseExtinction', responseExtinctionSeconds, ...
            'secondsToIncludeOffset', 0/1000 ...
     );
end

function rParams = updateEyeMovementParams(rParams, userParams)
      rParams.temporalParams  = modifyStructParams(rParams.temporalParams, ...
          'emPathType', userParams.emPathType ...
      );
end


% Method to update the optical image params
function rParams = updateOpticalImageParams(rParams, userParams, fieldOfVieDegs)
    rParams.oiParams = modifyStructParams(rParams.oiParams, ...
        	'blur', userParams.blur, ...
            'pupilDiamMm', userParams.pupilDiamMm, ...
            'wavefrontSpatialSamples', userParams.wavefrontSpatialSamples, ...
            'opticsModel', userParams.opticsModel,...
            'fieldOfViewDegs', max([userParams.minimumOpticalImagefieldOfViewDegs fieldOfVieDegs]) ...
        );
end

% Method to update the mosaic params
function rParams = updateMosaicParams(rParams, userParams, mosaicFOVdegs)
    if (isfield(rParams.mosaicParams, 'realisticSconeSubmosaic'))
        error('The realisticSconeSubmosaic has been removed. Pass directly sConeMinDistanceFactor and sConeFreeRadiusMicrons\n');
    end
        
    if (isempty(userParams.resamplingFactor))
        if (mosaicFOVdegs < 0.3)
            resamplingFactor = 13;
        elseif (mosaicFOVdegs < 0.5)
            resamplingFactor = 11;
        elseif (mosaicFOVdegs < 1.0)
            resamplingFactor = 9;
         elseif (mosaicFOVdegs < 2.0)
            resamplingFactor = 7;
        else
            resamplingFactor = 5;
        end
    else
        resamplingFactor = userParams.resamplingFactor;
    end
    
    if (isempty(userParams.maxGridAdjustmentIterations))
        if (mosaicFOVdegs < 0.5)
            maxGridAdjustmentIterations = 8000;
        elseif (mosaicFOVdegs < 1.0)
            maxGridAdjustmentIterations = 7500;
        elseif (mosaicFOVdegs < 2.0)
            maxGridAdjustmentIterations = 7000;
        elseif (mosaicFOVdegs < 4.0)
            maxGridAdjustmentIterations = 6500;
        elseif (mosaicFOVdegs < 6.0)
            maxGridAdjustmentIterations = 6000;
        elseif (mosaicFOVdegs < 8.0)
            maxGridAdjustmentIterations = 5500;
        elseif (mosaicFOVdegs < 10.0)
            maxGridAdjustmentIterations = 5000;
        else
            maxGridAdjustmentIterations = 4000;
        end
    else
        maxGridAdjustmentIterations = userParams.maxGridAdjustmentIterations;
    end
    
    rParams.mosaicParams = modifyStructParams(rParams.mosaicParams, ...
            'conePacking', userParams.conePacking, ... 
            'resamplingFactor', resamplingFactor, ...
            'fieldOfViewDegs', mosaicFOVdegs, ... 
            'innerSegmentSizeMicrons',userParams.innerSegmentSizeMicrons, ...
            'coneSpacingMicrons', userParams.coneSpacingMicrons, ...
            'apertureBlur',userParams.apertureBlur, ...
            'mosaicRotationDegs',userParams.mosaicRotationDegs,...
            'coneDarkNoiseRate',userParams.coneDarkNoiseRate,...
            'LMSRatio',userParams.LMSRatio,...
            'integrationTimeInSeconds', userParams.integrationTime, ...
            'isomerizationNoise', 'random',...                  % select from {'random', 'frozen', 'none'}
            'osNoise', 'random', ...                            % select from {'random', 'frozen', 'none'}
            'osModel', 'Linear');
        
    if strcmp(userParams.conePacking, 'hex')
        rParams.mosaicParams = modifyStructParams(rParams.mosaicParams, ...
            'eccBasedConeQuantalEfficiency', userParams.eccBasedConeQuantalEfficiency, ...
            'sConeMinDistanceFactor', userParams.sConeMinDistanceFactor, ...
            'sConeFreeRadiusMicrons', userParams.sConeFreeRadiusMicrons, ...
            'latticeAdjustmentPositionalToleranceF', userParams.latticeAdjustmentPositionalToleranceF, ...
            'latticeAdjustmentDelaunayToleranceF', userParams.latticeAdjustmentDelaunayToleranceF, ...
            'maxGridAdjustmentIterations', maxGridAdjustmentIterations, ...
            'marginF', userParams.marginF);
    end

    % Make sure mosaic noise parameters match the frozen noise flag passed in
    if (userParams.freezeNoise)
        if (strcmp(rParams.mosaicParams.isomerizationNoise, 'random'))
            rParams.mosaicParams.isomerizationNoise = 'frozen';
        end
        if (strcmp(rParams.mosaicParams.osNoise, 'random'))
            rParams.mosaicParams.osNoise = 'frozen';
        end
    end
        
end

% Method to generate test direction params
function testDirectionParams = generateTestDirectionParams(userParams)
    % Parameters that define the LM instances we'll generate here
    %
    % Use default LMPlane.
    testDirectionParams = instanceParamsGenerate;
    
    switch (userParams.coneContrastDirection)
        case 'L+M'
            testDirectionParams = modifyStructParams(testDirectionParams, ... 
                'instanceType', 'LMPlane', ...
                'startAngle', 45, ...
                'deltaAngle', 90, ...
                'nAngles', 1);
    
        case 'L+M+S'
            testDirectionParams = modifyStructParams(testDirectionParams, ... 
                'instanceType', 'LMSPlane');
            % Settings for L+M+S
            testDirectionParams.startAzimuthAngle = 45;
            testDirectionParams.deltaAzimuthAngle = 90;
            testDirectionParams.nAzimuthAngles = 1;
            testDirectionParams.startElevationAngle = 35.26;
            testDirectionParams.deltaElevationAngle = 90;
            testDirectionParams.nElevationAngles = 1;
        otherwise
            error('Cone contrast direction ''%s'' is not implemented', userParams.coneContrastDirection)
    end
    

    testDirectionParams = modifyStructParams(testDirectionParams, ...    
        'trialsNum', userParams.nTrainingSamples, ...
        'nContrastsPerDirection', userParams.nContrastsPerDirection, ...
        'lowContrast', userParams.lowContrast, ...
        'highContrast', userParams.highContrast, ...
        'contrastScale', userParams.contrastScale ...
        ); 
end

function thresholdParams = generateThresholdParams(userParams)
    % Parameters related to how we find thresholds from responses
    thresholdParams = thresholdParamsGenerate;
    thresholdParams = modifyStructParams(thresholdParams, ...
        'criterionFraction', userParams.thresholdCriterionFraction, ...
        'method', userParams.performanceClassifier, ...
        'STANDARDIZE', false, ...                   % Standardize data for PCA
        'standardizeSVMpredictors', false, ...       % Standardize data for SVM
        'useRBFKernel', userParams.useRBFSVMKernel, ...
        'PCAComponents', userParams.thresholdPCA, ...
        'signalSource', userParams.performanceSignal ...
        );
    
    if (strcmp(thresholdParams.method, 'svm')) || ...
        (strcmp(thresholdParams.method, 'svmV1FilterBank')) || ...
        (strcmp(thresholdParams.method, 'svmV1FilterEnsemble')) || ...
        (strcmp(thresholdParams.method, 'svmGaussianRF'))
        thresholdParams = modifyStructParams(thresholdParams, ...
            'spatialPoolingKernelParams', userParams.spatialPoolingKernelParams ...
            );  
    end

    if (isempty(userParams.performanceTrialsUsed))
        thresholdParams.trialsUsed = userParams.nTrainingSamples;
    else
        thresholdParams.trialsUsed = userParams.performanceTrialsUsed;
    end 
end



