function c_DavilaGeislerReplicate(varargin)
% c_DavilaGeislerReplicate(varargin)
%
% Compute thresholds to replicate spatial summation calculations of Davila and Geisler, more or less.
%
% This looks at thresholds as a function of spot size.  Our stimuli are
% monochromatic rather than monitor based, but to first order that should not make much difference
%
% Key/value pairs
%   'nTrainingSamples' - value (default 500).  Number of training samples to cycle through.
%   'spotDiametersMinutes' - vector (default [0.5 1 5 10 20 40]). Diameters of spots to be investigated.
%   'backgroundSizeDegs' - value (default 2.1). Size of square background
%   'wavelength' - value (default 550). Wavelength to use in calculations
%   'luminances' - vector (default [10]).  Background luminances in cd/m2 to be investigated.
%   'blur' - true/false (default true). Incorporate lens blur.
%   'imagePixels' - value (default 400).  Size of image pixel array
%   'computeResponses' - true/false (default true).  Compute responses.
%   'findPerformance' - true/false (default true).  Find performance.
%   'fitPsychometric' - true/false (default true).  Fit psychometric functions.
%   'plotPsychometric' - true/false (default true).  Plot psychometric functions.
%   'plotCSF' - true/false (default true).  Plot results.

%% Parse input
p = inputParser;
p.addParameter('nTrainingSamples',500,@isnumeric);
p.addParameter('spotDiametersMinutes',[0.5 1 5 10 20 40],@isnumeric);
p.addParameter('backgroundSizeDegs',[2.1],@isnumeric);
p.addParameter('wavelength',550,@isnumeric);
p.addParameter('luminances',[10],@isnumeric);
p.addParameter('blur',true,@islogical);
p.addParameter('imagePixels',400,@isnumeric);
p.addParameter('computeResponses',true,@islogical);
p.addParameter('findPerformance',true,@islogical);
p.addParameter('fitPsychometric',true,@islogical);
p.addParameter('plotPsychometric',true,@islogical);
p.addParameter('plotCSF',true,@islogical);
p.parse(varargin{:});

%% Get the parameters we need
%
% Start with default
rParams = responseParamsGenerate('spatialType','spot','backgroundType','AO','modulationType','AO');

%% Loop over spatial frequency
for ll = 1:length(p.Results.luminances)
    for cc = 1:length(p.Results.spotDiametersMinutes)
        
        % Get stimulus parameters correct
        %
        % Spatial 
        rParams.spatialParams.spotSizeDegs = p.Results.spotDiametersMinutes(cc)/60;
        rParams.spatialParams.backgroundSizeDegs = p.Results.backgroundSizeDegs;
        spatialParams.fieldOfViewDegs = 1.1*p.Results.backgroundSizeDegs;
        spatialParams.row = p.Results.imagePixels;
        spatialParams.col = p.Results.imagePixels;
        spatialParams.viewingDistance = 7.16;
        
        % Blur
        rParams.oiParams.blur = p.Results.blur;
        
        % Keep mosaic size in lock step with stimulus.  This is also forced before
        % the mosaic is created, but we need it here so that filenames are
        % consistent.  It is possible that we should not have a separate mosaic
        % size field, and just alwasy force it to match the scene.
        rParams.mosaicParams.fieldOfViewDegs = rParams.spatialParams.fieldOfViewDegs;
        
        % Set background luminance
        %
        % We start with a base luminance that we know is about mid-gray on the
        % monitor we specify.  To change luminance, we specify a scale factor.
        % This is eventually applied both to the background luminance and to the
        % monitor channel spectra, so that we don't get unintersting out of gamut errors.
        rParams.backgroundParams.backgroundWavelengthsNm = [p.Results.wavelength];
        rParams.backgroundParams.backgroundCornealPowerUW = [1];
        
        % Spot color parameters
        colorModulationParams.startWl = p.Results.wavelength;
        colorModulationParams.endWl = p.Results.wavelength;
        colorModulationParams.deltaWl = 10;
        colorModulationParams.spotWavelengthNm = p.Results.wavelength0;
        colorModulationParams.spotCornealPowerUW = 20;
        colorModulationParams.contrast = 1;
        
        % Pupil size.  They used a 3mm artificial pupil
        oiParams.pupilDiamMm = 3;
        
        % Set duration equal to sampling interval to do just one frame.
        %
        % Their intervals were 100 msec each.
        rParams.temporalParams.simulationTimeStepSecs = 100/1000;
        rParams.temporalParams.stimulusDurationInSeconds = rParams.temporalParams.simulationTimeStepSecs;
        rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.simulationTimeStepSecs;
        rParams.temporalParams.secondsToInclude = rParams.temporalParams.simulationTimeStepSecs;
        
        % Their main calculation was without eye movements
        rParams.temporalParams.eyesDoNotMove = true;
        
        % Set up mosaic parameters for just one stimulus time step
        rParams.mosaicParams.timeStepInSeconds = rParams.temporalParams.simulationTimeStepSecs;
        rParams.mosaicParams.integrationTimeInSeconds = rParams.mosaicParams.timeStepInSeconds;
        rParams.mosaicParams.isomerizationNoise = true;
        rParams.mosaicParams.osNoise = true;
        rParams.mosaicParams.osModel = 'Linear';
        
        % Parameters that define the LM instances we'll generate here
        %
        % Use default LMPlane.
        testDirectionParams = instanceParamsGenerate;
        testDirectionParams.startAngle = 45;
        testDirectionParams.deltaAngle = 90;
        testDirectionParams.nAngles = 1;
        
        % Number of contrasts to run in each color direction
        testDirectionParams.nContrastsPerDirection = 20;
        testDirectionParams.lowContrast = 0.0001;
        testDirectionParams.highContrast = 0.1;
        testDirectionParams.contrastScale = 'log';    % choose between 'linear' and 'log'
        
        % Parameters related to how we find thresholds from responses
        %
        % Use default
        thresholdParams = thresholdParamsGenerate;
        
        % Set number of trials
        testDirectionParams.trialsNum = p.Results.nTrainingSamples;
        
        %% Compute response instances
        if (p.Results.computeResponses)
            t_coneCurrentEyeMovementsResponseInstances('rParams',rParams,'testDirectionParams',testDirectionParams,'compute',true,'visualizeResponses',false);
        end
        
        %% Find performance, template max likeli
        thresholdParams.method = 'mlpt';
        if (p.Results.findPerformance)
            t_colorDetectFindPerformance('rParams',rParams,'testDirectionParams',testDirectionParams,'thresholdParams',thresholdParams,'compute',true,'plotSvmBoundary',false,'plotPsychometric',false);
        end
        
        %% Fit psychometric functions
        if (p.Results.fitPsychometric)
            davilaGeislerReplicate.spotDiametersMinutes(ll,cc) = p.Results.spotDiametersMinutes(cc);
            thresholdParams.method = 'mlpt';
            davilaGeislerReplicate.mlptThresholds(ll,cc) = t_plotDetectThresholdsOnLMPlane('rParams',rParams,'instanceParams',testDirectionParams,'thresholdParams',thresholdParams, ...
                'plotPsychometric',p.Results.plotPsychometric,'plotEllipse',false);
            close all;
        end
    end
end

%% Write out the data
%
% Set key spatial params to 0 to define a summary directory name
if (p.Results.fitPsychometric)
    fprintf('Writing performance data ... ');
    nameParams = rParams.spatialParams;
    nameParams.spotDiametersMinutes = 0;
    nameParams.fieldOfViewDegs = 0;
    nameParams.gaussianFWHMDegs = 0;
    paramsList = {nameParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    rwObject.write('davilaGeislerReplicate',davilaGeislerReplicate,paramsList,writeProgram);
    fprintf('done\n');
end

%% Make a plot of estimated threshold versus training set size
%
% The way the plot is coded counts on the test contrasts never changing
% across the conditions, which we could explicitly check for here.
if (p.Results.plotCSF)
    fprintf('Reading performance data ...');
    nameParams = rParams.spatialParams;
    nameParams.spotDiametersMinutes = 0;
    nameParams.fieldOfViewDegs = 0;
    nameParams.gaussianFWHMDegs = 0;
    paramsList = {nameParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    davilaGeislerReplicate = rwObject.read('davilaGeislerReplicate',paramsList,writeProgram);
    fprintf('done\n');
    
    hFig = figure; clf; hold on
    fontBump = 4;
    markerBump = -4;
    set(gcf,'Position',[100 100 450 650]);
    set(gca,'FontSize', rParams.plotParams.axisFontSize+fontBump);
    theColors = ['r' 'g' 'b'];
    legendStr = cell(length(p.Results.luminances),1);
    for ll = 1:length(p.Results.luminances)
        theColorIndex = rem(ll,length(theColors)) + 1;
        plot(davilaGeislerReplicate.spotDiametersMinutes(ll,:),1./[davilaGeislerReplicate.mlptThresholds(ll,:).thresholdContrasts]*davilaGeislerReplicate.mlptThresholds(1).testConeContrasts(1), ...
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
        rwObject.write('davilaGeislerReplicateWithBlur',hFig,paramsList,writeProgram,'Type','figure');
    else
        title(sprintf('Computational Observer CSF - no blur',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
        rwObject.write('davilaGeislerReplicateNoBlur',hFig,paramsList,writeProgram,'Type','figure');
    end
end

