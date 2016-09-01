function c_BanksEtAlReplicate(varargin)
% c_BanksEtAlReplicate(varargin)
%
% Compute thresholds to replicate Banks et al, 1987, more or less.
%
% This looks at L+M detection thrsholds, which seems close enough for right now to
% the isochromatic thresholds studied by Banks et al.
%
% Key/value pairs
%   'nTrainingSamples' - value (default 500).  Number of training samples to cycle through.
%   'cyclesPerDegree' - vector (default [3 5 10 20 40]). Spatial frequencoes of grating to be investigated.
%   'luminances' - vector (default [3.4 34 340]).  Luminances in cd/m2 to be investigated.
%   'blur' - true/false (default true). Incorporate lens blur.
%   'computeResponses' - true/false (default true).  Compute responses.
%   'fitPsychometric' - true/false (default true).  Fit psychometric functions.
%   'plotPsychometric' - true/false (default true).  Plot psychometric functions.
%   'plotCSF' - true/false (default true).  Plot results.

%% Parse input
p = inputParser;
p.addParameter('nTrainingSamples',500,@isnumeric);
p.addParameter('cyclesPerDegree',[3 5 10 20 40 50],@isnumeric);
p.addParameter('luminances',[3.4 34 340],@isnumeric);
p.addParameter('blur',true,@islogical);
p.addParameter('computeResponses',true,@islogical);
p.addParameter('fitPsychometric',true,@islogical);
p.addParameter('plotPsychometric',true,@islogical);
p.addParameter('plotCSF',true,@islogical);
p.parse(varargin{:});

%% Get the parameters we need
%
% Start with default
rParams = colorGaborResponseParamsGenerate;

%% Loop over spatial frequency
for ll = 1:length(p.Results.luminances)
    for cc = 1:length(p.Results.cyclesPerDegree)
        
        % Get stimulus parameters correct
        %
        % The stimulus was half-cosine windowed to contain 7.5 cycles.  We set
        % our half-cosine window to match that and also make the field of view
        % just a tad bigger.
        rParams.gaborParams.windowType = 'halfcos';
        rParams.gaborParams.cyclesPerDegree = p.Results.cyclesPerDegree(cc);
        rParams.gaborParams.gaussianFWHMDegs = 3.75*(1/rParams.gaborParams.cyclesPerDegree);
        rParams.gaborParams.fieldOfViewDegs = 2.1*rParams.gaborParams.gaussianFWHMDegs;
        
        % Blur
        rParams.oiParams.blur = p.Results.blur;
        
        % Keep mosaic size in lock step with stimulus.  This is also forced before
        % the mosaic is created, but we need it here so that filenames are
        % consistent.  It is possible that we should not have a separate mosaic
        % size field, and just alwasy force it to match the scene.
        rParams.mosaicParams.fieldOfViewDegs = rParams.gaborParams.fieldOfViewDegs;
        
        % Set background luminance
        %
        % We start with a base luminance that we know is about mid-gray on the
        % monitor we specify.  To change luminance, we specify a scale factor.
        % This is eventually applied both to the background luminance and to the
        % monitor channel spectra, so that we don't get unintersting out of gamut errors.
        baseLum = 50;
        theLum = p.Results.luminances(ll);
        rParams.backgroundParams.backgroundxyY = [0.33 0.33 baseLum]';
        rParams.backgroundParams.monitorFile = 'CRT-MODEL';
        rParams.backgroundParams.leakageLum = 1.0;
        rParams.backgroundParams.lumFactor = theLum/baseLum;
        
        % Pupil size.  They used a 2mm artificial pupil
        oiParams.pupilDiamMm = 2;
        
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
        testDirectionParams = LMPlaneInstanceParamsGenerate;
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
            t_colorGaborConeCurrentEyeMovementsResponseInstances('rParams',rParams,'testDirectionParams',testDirectionParams,'compute',true,'visualizeResponses',false);
        end
        
        %% Find performance, template max likeli
        thresholdParams.method = 'mlpt';
        t_colorGaborDetectFindPerformance('rParams',rParams,'testDirectionParams',testDirectionParams,'thresholdParams',thresholdParams,'compute',true,'plotSvmBoundary',false,'plotPsychometric',false);
        
        %% Fit psychometric functions
        if (p.Results.fitPsychometric)
            banksEtAlReplicate.cyclesPerDegree(ll,cc) = p.Results.cyclesPerDegree(cc);
            thresholdParams.method = 'mlpt';
            banksEtAlReplicate.mlptThresholds(ll,cc) = t_plotGaborDetectThresholdsOnLMPlane('rParams',rParams,'LMPlaneInstanceParams',testDirectionParams,'thresholdParams',thresholdParams, ...
                'plotPsychometric',p.Results.plotPsychometric,'plotEllipse',false);
            close all;
        end
    end
end

%% Write out the data
%
% Set cyclesPerDegree to 0 to define a summary directory name
if (p.Results.fitPsychometric)
    fprintf('Writing performance data ... ');
    rParams.gaborParams.cyclesPerDegree = 0;
    paramsList = {rParams.gaborParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    rwObject.write('banksEtAlReplicate',banksEtAlReplicate,paramsList,writeProgram);
    fprintf('done\n');
end

%% Make a plot of estimated threshold versus training set size
%
% The way the plot is coded counts on the test contrasts never changing
% across the conditions, which we could explicitly check for here.
if (p.Results.plotCSF)
    fprintf('Reading performance data ...');
    rParams.gaborParams.cyclesPerDegree = 0;
    paramsList = {rParams.gaborParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams};
    rwObject = IBIOColorDetectReadWriteBasic;
    writeProgram = mfilename;
    banksEtAlReplicate = rwObject.read('banksEtAlReplicate',paramsList,writeProgram);
    fprintf('done\n');
    
    hFig = figure; clf; hold on
    fontBump = 4;
    markerBump = -4;
    set(gcf,'Position',[100 100 450 650]);
    set(gca,'FontSize', rParams.plotParams.axisFontSize+fontBump);
    theColors = ['r' 'g' 'b'];
    for ll = 1:length(p.Results.luminances)
        theColorIndex = rem(ll,length(theColors)) + 1;
        plot(banksEtAlReplicate.cyclesPerDegree,1./[banksEtAlReplicate.mlptThresholds(ll,:).thresholdContrasts]*banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1), ...
            [theColors(theColorIndex) 'o'],'MarkerSize',rParams.plotParams.markerSize+markerBump,'MarkerFaceColor',theColors(theColorIndex));
        plot(banksEtAlReplicate.cyclesPerDegree,1./[banksEtAlReplicate.mlptThresholds(ll,:).thresholdContrasts]*banksEtAlReplicate.mlptThresholds(1).testConeContrasts(1), ...
            theColors(theColorIndex),'LineWidth',rParams.plotParams.lineWidth);
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
        title(sprintf('Computational Observer CSF - w/ Blur',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
    else
        title(sprintf('Computational Observer CSF - no/ Blur',rParams.mosaicParams.fieldOfViewDegs'),'FontSize',rParams.plotParams.titleFontSize+fontBump);
    end
    rwObject.write('banksEtAlReplicate',hFig,paramsList,writeProgram,'Type','figure');
end

