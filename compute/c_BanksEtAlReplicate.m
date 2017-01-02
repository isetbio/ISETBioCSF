function [validationData, extraData] = c_BanksEtAlReplicate(varargin)
% [validationData, extraData] = c_BanksEtAlReplicate(varargin)
%
% Compute thresholds to replicate Banks et al, 1987, more or less.  That
% paper builds on methods reported in Geisler, 1984.
%
% This looks at L+M detection thrsholds, which seems close enough for right now to
% the isochromatic thresholds studied by Banks et al.
%
% Key/value pairs
%   'useScratchTopLevelDirName'- true/false (default false). 
%      When true, the top level output directory is [scratch]. 
%      When false, it is the name of this script.
%   'nTrainingSamples' - value (default 500).  Number of training samples to cycle through.
%   'cyclesPerDegree' - vector (default [3 5 10 20 40]). Spatial frequencoes of grating to be investigated.
%   'luminances' - vector (default [3.4 34 340]).  Luminances in cd/m2 to be investigated.
%   'pupilDiamMm' - value (default 2).  Pupil diameter in mm.
%   'blur' - true/false (default true). Incorporate lens blur.
%   'innerSegmentSizeMicrons' - Diameter of the cone light-collecting area, in microns 
%       Default: sizeForSquareApertureFromDiameterForCircularAperture(3.0), where 3 microns = 6 min arc for 300 mirons/degree in the human retina.
%   'apertureBlur' - Blur by cone aperture? true/false (default false).
%   'conePacking'   - how cones are packed spatially. 
%       Choose from : 'rect', for a rectangular mosaic
%                     'hex', for a hex mosaic with an eccentricity-varying cone spacing
%                     'hexReg' for a hex mosaic witha regular cone spacing
%   'imagePixels' - value (default 400).  Size of image pixel array
%   'nContrastsPerDirection' - value (default 20). Number of contrasts.
%   'lowContrast' - value (default 0.0001). Low contrast.
%   'highContrast' - value (default 0.1). High contrast.
%   'contrastScale' - 'log'/'linear' (default 'log'). Contrast scale.
%   'computeResponses' - true/false (default true).  Compute responses.
%   'findPerformance' - true/false (default true).  Find performance.
%   'fitPsychometric' - true/false (default true).  Fit psychometric functions.
%   'thresholdCriterionFraction' value (default 0.75). Criterion corrrect for threshold.
%   'generatePlots' - true/false (default true).  No plots are generated unless this is true.
%   'visualizedResponseNormalization' - string (default 'submosaicBasedZscore'). How to normalize visualized responses
%        Available options: 'submosaicBasedZscore', 'LMSabsoluteResponseBased', 'LMabsoluteResponseBased', 'MabsoluteResponseBased'
%   'plotPsychometric' - true/false (default true).  Plot psychometric functions.
%   'plotCSF' - true/false (default true).  Plot results.
%   'freezeNoise' - true/false (default true). Freeze noise so calculations reproduce.

%% Parse input
p = inputParser;
p.addParameter('useScratchTopLevelDirName', false, @islogical);
p.addParameter('nTrainingSamples',500,@isnumeric);
p.addParameter('cyclesPerDegree',[3 5 10 20 40 50],@isnumeric);
p.addParameter('luminances',[3.4 34 340],@isnumeric);
p.addParameter('pupilDiamMm',2,@isnumeric);
p.addParameter('blur',true,@islogical);
p.addParameter('innerSegmentSizeMicrons',sizeForSquareApertureFromDiameterForCircularAperture(3.0), @isnumeric);   % 3 microns = 0.6 min arc for 300 microns/deg in human retina
p.addParameter('apertureBlur', false, @islogical);
p.addParameter('coneSpacingMicrons', 3.0, @isnumeric);
p.addParameter('conePacking', 'hexReg');                 
p.addParameter('imagePixels',400,@isnumeric);
p.addParameter('nContrastsPerDirection',20,@isnumeric);
p.addParameter('lowContrast',0.0001,@isnumeric);
p.addParameter('highContrast',0.1,@isnumeric);
p.addParameter('contrastScale','log',@ischar);
p.addParameter('computeResponses',true,@islogical);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.addParameter('findPerformance',true,@islogical);
p.addParameter('fitPsychometric',true,@islogical);
p.addParameter('thresholdCriterionFraction',0.75,@isnumeric);
p.addParameter('generatePlots',true,@islogical);
p.addParameter('plotPsychometric',true,@islogical);
p.addParameter('plotCSF',true,@islogical);
p.addParameter('freezeNoise',true,@islogical);
p.parse(varargin{:});

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
        
        % Get stimulus parameters correct
        %
        % The stimulus was half-cosine windowed to contain 7.5 cycles.  We set
        % our half-cosine window to match that and also make the field of view
        % just a tad bigger.
        cyclesPerDegree = p.Results.cyclesPerDegree(cc);
        gaussianFWHMDegs = 3.75*(1/cyclesPerDegree);
        fieldOfViewDegs = 2.1*gaussianFWHMDegs;
        rParams.spatialParams = modifyStructParams(rParams.spatialParams, ...
            'windowType', 'halfcos', ...
            'cyclesPerDegree', cyclesPerDegree, ...
            'gaussianFWHMDegs', gaussianFWHMDegs, ...
            'fieldOfViewDegs', fieldOfViewDegs, ...
            'row', p.Results.imagePixels, ...
            'col', p.Results.imagePixels);
        
        % Blur
        %
        % Banks et al. used a 2mm artificial pupil
        rParams.oiParams = modifyStructParams(rParams.oiParams, ...
        	'blur', p.Results.blur, ...
            'pupilDiamMm', p.Results.pupilDiamMm ...    
        );
              
        % Set background luminance
        %
        % We start with a base luminance that we know is about mid-gray on the
        % monitor we specify.  To change luminance, we specify a scale factor.
        % This is eventually applied both to the background luminance and to the
        % monitor channel spectra, so that we don't get unintersting out of gamut errors.
        baseLum = 50;
        theLum = p.Results.luminances(ll);
        rParams.backgroundParams = modifyStructParams(rParams.backgroundParams, ...
        	'backgroundxyY', [0.33 0.33 baseLum]',...
        	'monitorFile', 'CRT-MODEL', ...
        	'leakageLum', 1.0, ...
        	'lumFactor', theLum/baseLum);
        
        % Their stimulus intervals were 100 msec each.
        %
        % Equate stimulusSamplingIntervalInSeconds to stimulusDurationInSeconds to generate 1 time point only.
        % Their main calculation was without eye movements
        stimulusDurationInSeconds = 100/1000;
        rParams.temporalParams = modifyStructParams(rParams.temporalParams, ...
            'stimulusDurationInSeconds', stimulusDurationInSeconds, ...
            'stimulusSamplingIntervalInSeconds',  stimulusDurationInSeconds, ... 
            'secondsToInclude', stimulusDurationInSeconds, ...
            'emPathType', 'none' ...       
        );
        
        % Set up mosaic parameters. Here we integrate for the entire stimulus duration (100/1000)
        %
        % Keep mosaic size in lock step with stimulus
        rParams.mosaicParams = modifyStructParams(rParams.mosaicParams, ...
            'fieldOfViewDegs', rParams.spatialParams.fieldOfViewDegs, ...  
            'innerSegmentSizeMicrons',p.Results.innerSegmentSizeMicrons, ...
            'apertureBlur',p.Results.apertureBlur, ...
            'coneSpacingMicrons', p.Results.coneSpacingMicrons, ...
            'conePacking', p.Results.conePacking, ...
        	'integrationTimeInSeconds', rParams.temporalParams.stimulusDurationInSeconds, ...
        	'isomerizationNoise', 'random',...              % select from {'random', 'frozen', 'none'}
        	'osNoise', 'random', ...                        % select from {'random', 'frozen', 'none'}
        	'osModel', 'Linear');
        
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
        % Use default
        thresholdParams = thresholdParamsGenerate;
        thresholdParams.criterionFraction = p.Results.thresholdCriterionFraction;
        
        %% Compute response instances
        desiredTrialsPerBlock = 50;
        trialBlocks = round(testDirectionParams/desiredTrialsPerBlock);
        if (p.Results.computeResponses)
           t_coneCurrentEyeMovementsResponseInstances(...
               'rParams',rParams,...
               'testDirectionParams',testDirectionParams,...
               'compute',true,...
               'visualizedResponseNormalization', p.Results.visualizedResponseNormalization, ...
               'generatePlots',p.Results.generatePlots,...
               'freezeNoise',p.Results.freezeNoise,...
               'trialBlocks',trialBlocks);
        end
        
        %% Find performance, template max likeli
        thresholdParams.method = 'mlpt';
        if (p.Results.findPerformance)
            t_colorDetectFindPerformance('rParams',rParams,'testDirectionParams',testDirectionParams,'thresholdParams',thresholdParams,'compute',true,'plotSvmBoundary',false,'plotPsychometric',false,'freezeNoise',p.Results.freezeNoise);
        end
        
        %% Fit psychometric functions
        if (p.Results.fitPsychometric)
            banksEtAlReplicate.cyclesPerDegree(ll,cc) = p.Results.cyclesPerDegree(cc);
            thresholdParams.method = 'mlpt';
            banksEtAlReplicate.mlptThresholds(ll,cc) = t_plotDetectThresholdsOnLMPlane('rParams',rParams,'instanceParams',testDirectionParams,'thresholdParams',thresholdParams, ...
                'plotPsychometric',p.Results.generatePlots & p.Results.plotPsychometric,'plotEllipse',false);
            %close all;
        end
    end
end

%% Write out the data
if (p.Results.fitPsychometric)
    fprintf('Writing performance data ... ');
    paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams,  rParams.backgroundParams, testDirectionParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    rwObject.write('banksEtAlReplicate',banksEtAlReplicate,paramsList,writeProgram);
    fprintf('done\n');
end

%% Get performance data
fprintf('Reading performance data ...');
paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams,  rParams.backgroundParams, testDirectionParams};
rwObject = IBIOColorDetectReadWriteBasic;
writeProgram = mfilename;
banksEtAlReplicate = rwObject.read('banksEtAlReplicate',paramsList,writeProgram);
fprintf('done\n');

%% Output validation data
if (nargout > 0)
    validationData.cyclesPerDegree = banksEtAlReplicate.cyclesPerDegree;
    validationData.mlptThresholds = banksEtAlReplicate.mlptThresholds;
    validationData.luminances = p.Results.luminances;
    extraData.paramsList = paramsList;
    extraData.p.Results = p.Results;
end

%% Make a plot of estimated threshold versus training set size
%
% The way the plot is coded counts on the test contrasts never changing
% across the conditions, which we could explicitly check for here.
if (p.Results.generatePlots & p.Results.plotCSF)  
    hFig = figure; clf; hold on
    fontBump = 4;
    markerBump = -4;
    set(gcf,'Position',[100 100 450 650]);
    set(gca,'FontSize', rParams.plotParams.axisFontSize+fontBump);
    theColors = ['r' 'g' 'b'];
    legendStr = cell(length(p.Results.luminances),1);
    for ll = 1:length(p.Results.luminances)
        theColorIndex = rem(ll,length(theColors)) + 1;
        plot(banksEtAlReplicate.cyclesPerDegree(ll,:),1./[banksEtAlReplicate.mlptThresholds(ll,:).thresholdContrasts]*banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1), ...
            [theColors(theColorIndex) 'o-'],'MarkerSize',rParams.plotParams.markerSize+markerBump,'MarkerFaceColor',theColors(theColorIndex),'LineWidth',rParams.plotParams.lineWidth);  
        legendStr{ll} = sprintf('%0.1f cd/m2',p.Results.luminances(ll));
    end
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    xlabel('Log10 Spatial Frequency (cpd)', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
    ylabel('Log10 Contrast Sensitivity', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
    xlim([1 100]); ylim([10 10000]);
    legend(legendStr,'Location','NorthEast','FontSize',rParams.plotParams.labelFontSize+fontBump);
    box off; grid on
    if (p.Results.blur)
        title(sprintf('Computational Observer CSF - w/ blur',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
        rwObject.write('banksEtAlReplicateWithBlur',hFig,paramsList,writeProgram,'Type','figure');
    else
        title(sprintf('Computational Observer CSF - no blur',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
        rwObject.write('banksEtAlReplicateNoBlur',hFig,paramsList,writeProgram,'Type','figure');
    end
end

