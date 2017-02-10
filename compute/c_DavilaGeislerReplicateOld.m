function c_DavilaGeislerReplicate(varargin)
% c_DavilaGeislerReplicate(varargin)
%
% Compute thresholds to replicate spatial summation calculations of Davila and Geisler, more or less.
%
% This looks at thresholds as a function of spot size.  Our stimuli are
% monochromatic rather than monitor based, but to first order that should not make much difference
%
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
%   'imagePixels' - value (default 400).  Size of image pixel array
%   'computeResponses' - true/false (default true).  Compute responses.
%   'minContrast' - value (default 1e-6). Minimum contrast to use
%   'maxContrast' - value (default 1e-2). Max contrast to use
%   'nContrasts' - value (default 30). Number of contrasts to use
%   'findPerformance' - true/false (default true).  Find performance.
%   'fitPsychometric' - true/false (default true).  Fit psychometric functions.
%   'generatePlots' - true/false (default true).  Generate plots?  Other
%     plot options only have an effect if this is true.
%   'plotPsychometric' - true/false (default true).  Plot psychometric functions.
%   'plotSpatialSummation' - true/false (default true).  Plot results.

%% Parse input
p = inputParser;
p.addParameter('useScratchTopLevelDirName', false, @islogical);
p.addParameter('nTrainingSamples',500,@isnumeric);
p.addParameter('spotDiametersMinutes',[0.5 0.75 1 1.5 2 2.5 5 10 20 40],@isnumeric);
p.addParameter('backgroundSizeDegs',[2.1],@isnumeric);
p.addParameter('wavelength',550,@isnumeric);
p.addParameter('luminances',[10],@isnumeric);
p.addParameter('pupilDiamMm',3,@isnumeric);
p.addParameter('durationMs',100,@isnumeric);
p.addParameter('blur',true,@islogical);
p.addParameter('imagePixels',400,@isnumeric);
p.addParameter('computeResponses',true,@islogical);
p.addParameter('minContrast',1e-6,@isnumeric);
p.addParameter('maxContrast',1e-2,@isnumeric);
p.addParameter('nContrasts',30,@isnumeric);
p.addParameter('findPerformance',true,@islogical);
p.addParameter('fitPsychometric',true,@islogical);
p.addParameter('generatePlots',true,@islogical);
p.addParameter('plotPsychometric',true,@islogical);
p.addParameter('plotSpatialSummation',true,@islogical);
p.parse(varargin{:});

%% Get the parameters we need
%
% Start with default
rParams = responseParamsGenerate('spatialType','spot','backgroundType','AO','modulationType','AO');

%% Set the  topLevelDir name
if (~p.Results.useScratchTopLevelDirName)
    rParams.topLevelDirParams.name = mfilename;
end

%% Loop over spatial frequency
for ll = 1:length(p.Results.luminances)
    for dd = 1:length(p.Results.spotDiametersMinutes)
        
        % Get stimulus parameters correct
        %
        % Spatial 
        rParams.spatialParams.spotSizeDegs = p.Results.spotDiametersMinutes(dd)/60;
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
        
        % Pupil size.  They used a 3mm artificial pupil
        rParams.oiParams.pupilDiamMm = p.Results.pupilDiamMm;
        
        % Scale radiance to produce desired background levels
        deltaWl = 10;
        S = [p.Results.wavelength deltaWl 2];
        wls = SToWls(S);
        load T_xyz1931
        T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
        radiancePerUW = AOMonochromaticCornealPowerToRadiance(wls,rParams.backgroundParams.backgroundWavelengthsNm,1,rParams.oiParams.pupilDiamMm,rParams.spatialParams.backgroundSizeDegs^2);
        xyzPerUW = T_xyz*radiancePerUW;
        desiredLumCdM2 = p.Results.luminances(ll);
        rParams.backgroundParams.backgroundCornealPowerUW = desiredLumCdM2/xyzPerUW(2);
        
        % Scale spot parameters into what we think a reasonable range is.
        % In Davila & Geisler, they give thresholds in what they call
        % threshold energy, which is spotLuminance*durationSecs*areaMinutes2.
        %
        % Just to get things into a reasonable scale, we'll find a
        % luminance corresponding to a threshold energy of 10 for our
        % smallest spot.
        maxThresholdEnergy = 10;
        minSpotAreaMin2 = pi*(min(p.Results.spotDiametersMinutes)/2)^2;
        maxSpotLuminanceCdM2 = maxThresholdEnergy/((p.Results.durationMs/1000)*minSpotAreaMin2);
        rParams.colorModulationParams.startWl = p.Results.wavelength;
        rParams.colorModulationParams.endWl = p.Results.wavelength+deltaWl;
        rParams.colorModulationParams.deltaWl = deltaWl;
        rParams.colorModulationParams.spotWavelengthNm = p.Results.wavelength;
        rParams.colorModulationParams.spotCornealPowerUW = maxSpotLuminanceCdM2/xyzPerUW(2);
        rParams.colorModulationParams.contrast = 1;
     
        % Set duration equal to sampling interval to do just one frame.
        %
        % Their intervals were 100 msec each.
        rParams.temporalParams.simulationTimeStepSecs = p.Results.durationMs/1000;
        rParams.temporalParams.stimulusDurationInSeconds = rParams.temporalParams.simulationTimeStepSecs;
        rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.simulationTimeStepSecs;
        rParams.temporalParams.secondsToInclude = rParams.temporalParams.simulationTimeStepSecs;
        
        % Their main calculation was without eye movements
        rParams.temporalParams.eyesDoNotMove = true;
        
        % Set up mosaic parameters for just one stimulus time step
        rParams.mosaicParams.timeStepInSeconds = rParams.temporalParams.simulationTimeStepSecs;
        rParams.mosaicParams.integrationTimeInSeconds = rParams.mosaicParams.timeStepInSeconds;
        rParams.mosaicParams.isomerizationNoise = 'frozen';
        rParams.mosaicParams.osNoise = 'frozen';
        rParams.mosaicParams.osModel = 'Linear';
        
        % Parameters that define the contrasts we'll study here
        testDirectionParams = instanceParamsGenerate('instanceType','contrasts');
    
        % Number of contrasts to run in each color direction
        testDirectionParams.nContrastsPerDirection = p.Results.nContrasts;
        testDirectionParams.lowContrast = p.Results.minContrast;
        testDirectionParams.highContrast = p.Results.maxContrast;
        testDirectionParams.contrastScale = 'log';    % choose between 'linear' and 'log'
        testDirectionParams.trialsNum = p.Results.nTrainingSamples;

        % Parameters related to how we find thresholds from responses
        %
        % Use default
        thresholdParams = thresholdParamsGenerate;

        %% Compute response instances
        if (p.Results.computeResponses)
            t_coneCurrentEyeMovementsResponseInstances('rParams',rParams,'testDirectionParams',testDirectionParams,'compute',true,'generatePlots',false);
        end
        
        %% Find performance, template max likeli
        thresholdParams.method = 'mlpt';
        if (p.Results.findPerformance)
            t_colorDetectFindPerformance('rParams',rParams,'testDirectionParams',testDirectionParams,'thresholdParams',thresholdParams,'compute',true,'plotSvmBoundary',false,'plotPsychometric',false);
        end
        
        %% Fit psychometric functions
        if (p.Results.fitPsychometric)
            davilaGeislerReplicate.spotDiametersMinutes(ll,dd) = p.Results.spotDiametersMinutes(dd);
            thresholdParams.method = 'mlpt';
            davilaGeislerReplicate.mlptThresholds(ll,dd) = t_fitPsychometricFunctions('rParams',rParams,'instanceParams',testDirectionParams,'thresholdParams',thresholdParams, ...
                'generatePlots',p.Results.generatePlots && p.Results.plotPsychometric);
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
    nameParams.spotSizeDegs = 0;
    nameParams.backgroundSizeDegs = 0;
    nameParams.fieldOfViewDegs = 0;
    paramsList = {rParams.topLevelDirParams, nameParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    rwObject.write('davilaGeislerReplicate',davilaGeislerReplicate,paramsList,writeProgram);
    fprintf('done\n');
end

%% Make a plot of estimated threshold versus training set size
%
% The way the plot is coded counts on the test contrasts never changing
% across the conditions, which we could explicitly check for here.
if (p.Results.generatePlots && p.Results.plotSpatialSummation)
    fprintf('Reading performance data ...');
    nameParams = rParams.spatialParams;
    nameParams.spotSizeDegs = 0;
    nameParams.backgroundSizeDegs = 0;
    nameParams.fieldOfViewDegs = 0;
    paramsList = {rParams.topLevelDirParams, nameParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    davilaGeislerReplicate = rwObject.read('davilaGeislerReplicate',paramsList,writeProgram);
    fprintf('done\n');
    
    % Some numbers we need
    spotAreasMin2 = pi*((p.Results.spotDiametersMinutes/2).^2);
    maxThresholdEnergies = maxSpotLuminanceCdM2*rParams.temporalParams.stimulusDurationInSeconds*spotAreasMin2;
    
    hFig = figure; clf; hold on
    fontBump = 4;
    markerBump = -4;
    set(gcf,'Position',[100 100 450 650]);
    set(gca,'FontSize', rParams.plotParams.axisFontSize+fontBump);
    theColors = ['r' 'g' 'b'];
    legendStr = cell(length(p.Results.luminances),1);
    for ll = 1:length(p.Results.luminances)
        theColorIndex = rem(ll,length(theColors)) + 1;
        plot(spotAreasMin2,[davilaGeislerReplicate.mlptThresholds(ll,:).thresholdContrasts].*maxThresholdEnergies, ...
            [theColors(theColorIndex) 'o-'],'MarkerSize',rParams.plotParams.markerSize+markerBump,'MarkerFaceColor',theColors(theColorIndex),'LineWidth',rParams.plotParams.lineWidth);  
        legendStr{ll} = sprintf('%0.1f cd/m2',p.Results.luminances(ll));
    end
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    xlabel('Log10 Spot Area (square arc minutes)', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
    ylabel('Log10 Threshold Energy', 'FontSize' ,rParams.plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
    xlim([1e-1 1e4]); ylim([1e-3 1]);
    legend(legendStr,'Location','NorthWest','FontSize',rParams.plotParams.labelFontSize+fontBump);
    box off; grid on
    if (p.Results.blur)
        title(sprintf('Computational Observer CSF - w/ blur',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
        rwObject.write('davilaGeislerReplicateWithBlur',hFig,paramsList,writeProgram,'Type','figure');
    else
        title(sprintf('Computational Observer CSF - no blur',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
        rwObject.write('davilaGeislerReplicateNoBlur',hFig,paramsList,writeProgram,'Type','figure');
    end
end

