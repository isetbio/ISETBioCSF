function generateAllScenes(noiseInstances, stimulusSizeDegs, meanLuminanceCdPerM2, contrastLevels, resourcesDir)
   
    viewingDistance = 75/100;
    
    % Generate stimulus spatial modulations
    noiseNorm = nan;
    
    
    nContrasts = numel(contrastLevels);
    oriDegs = 90;
    
    
    [lowFrequencySpatialModulations, lowFrequencySpatialModulationsOrtho, lowFrequencyTemplate, lowFrequencyTemplateOrtho, spatialSupportDegs, noiseNorm] = ...
        generateStimulusSpatialModulation(stimulusSizeDegs, noiseNorm, 'low frequency', oriDegs, contrastLevels, noiseInstances);
    
    [highFrequencySpatialModulations, highFrequencySpatialModulationsOrtho, highFrequencyTemplate, highFrequencyTemplateOrtho, spatialSupportDegs, noiseNorm] = ...
        generateStimulusSpatialModulation(stimulusSizeDegs, noiseNorm, 'high frequency', oriDegs, contrastLevels, noiseInstances);
 
    % Visualize the spatial modulations
    visualizeSpatialModulations(lowFrequencySpatialModulations, lowFrequencySpatialModulationsOrtho, ...
                    highFrequencySpatialModulations, highFrequencySpatialModulationsOrtho);
                
    
    % Generate ISETBio display
    presentationDisplay = generateDisplay(viewingDistance, 4*meanLuminanceCdPerM2);
    
    % Generate the null scene (zero spatial modulation)
    nullModulation = 0.5+0*squeeze(lowFrequencySpatialModulations(1, 1,:,:,:));
    nullScene = generateScene(nullModulation, ...
                presentationDisplay, stimulusSizeDegs, meanLuminanceCdPerM2);
            
    % Generate ISETBio scenes from spatial modulations
    for theContrastLevel = 1:nContrasts
        for theInstance = 1:noiseInstances 
            
            lowFrequencyScenes{theContrastLevel, theInstance} = generateScene(...
                squeeze(lowFrequencySpatialModulations(theInstance, theContrastLevel,:,:,:)), ...
                presentationDisplay,stimulusSizeDegs, meanLuminanceCdPerM2);
            
            lowFrequencyScenesOrtho{theContrastLevel, theInstance} = generateScene(...
                squeeze(lowFrequencySpatialModulationsOrtho(theInstance, theContrastLevel,:,:,:)), ...
                presentationDisplay,stimulusSizeDegs, meanLuminanceCdPerM2);
             
            highFrequencyScenes{theContrastLevel, theInstance} = generateScene(...
                squeeze(highFrequencySpatialModulations(theInstance, theContrastLevel,:,:,:)), ...
                presentationDisplay,stimulusSizeDegs, meanLuminanceCdPerM2);
 
            highFrequencyScenesOrtho{theContrastLevel, theInstance} = generateScene(...
                squeeze(highFrequencySpatialModulationsOrtho(theInstance, theContrastLevel,:,:,:)), ...
                presentationDisplay,stimulusSizeDegs, meanLuminanceCdPerM2);
        end
    end
   
    fName = fullfile(resourcesDir, sprintf('scenes_luminance_%2.1f.mat', meanLuminanceCdPerM2));
    save(fName, 'nullScene', 'lowFrequencyScenes', 'highFrequencyScenes', ...
         'lowFrequencyScenesOrtho', 'highFrequencyScenesOrtho', ...
         'lowFrequencyTemplate', 'lowFrequencyTemplateOrtho', ...
         'highFrequencyTemplate', 'highFrequencyTemplateOrtho', ...
         'spatialSupportDegs', ...
         'contrastLevels', 'noiseInstances', '-v7.3');
end


function [stimulusSpatialModulation, stimulusSpatialModulationOrtho, template, templateOrtho, spatialSupportDegs, noiseNorm] = ...
    generateStimulusSpatialModulation(stimulusSizeDegs, noiseNorm, stimulusType, oriDegs, contrastLevels, noiseInstances)

    coneApertureMicrons = 2; micronsPerDegree = 300;
    stimulusPixelSizeArcMin = 0.75*coneApertureMicrons / micronsPerDegree * 60;
    stimulusWidthArcMin = stimulusSizeDegs * 60;
    
    % Grating params
    gratingParams.oriDegs = oriDegs;
    gratingParams.sigmaArcMin = stimulusSizeDegs/7*60;
    gratingParams.contrastLevels = contrastLevels;
    % Noise params
    noiseParams.steepness = 100;

    switch (stimulusType)
        case 'low frequency'
            gratingParams.sfCPD = 4;
            noiseParams.spectrumShape = 'highPassCornerFreq';
            noiseParams.cornerFrequencyCPD = 10;
        case 'high frequency'
            gratingParams.sfCPD = 11;
            noiseParams.spectrumShape = 'lowPassCornerFreq';
            noiseParams.cornerFrequencyCPD = 5;
        otherwise
            error('Unknown stimulus type: ''%s''.', stimulusType);
    end
    
    
    figNo = 1;
    [stimStruct, noiseNorm] = ...
        generateGratingInNoiseSpatialModulationPattern(...
        stimulusWidthArcMin, stimulusPixelSizeArcMin,  ...
        gratingParams, noiseParams, noiseNorm, noiseInstances, figNo);
    
    spatialSupportDegs = stimStruct.spatialSupportDegs;
    
    % Normalize to [0 .. 1]
    stimulusSpatialModulation(:,:,:,:,1) = (squeeze(stimStruct.image(:,:,:,:)) + 1)/2;
    stimulusSpatialModulation(:,:,:,:,2) = stimulusSpatialModulation(:,:,:,:,1);
    stimulusSpatialModulation(:,:,:,:,3) = stimulusSpatialModulation(:,:,:,:,1);
    
    stimulusSpatialModulationOrtho(:,:,:,:,1) = (squeeze(stimStruct.imageOrtho(:,:,:,:)) + 1)/2;
    stimulusSpatialModulationOrtho(:,:,:,:,2) = stimulusSpatialModulationOrtho(:,:,:,:,1);
    stimulusSpatialModulationOrtho(:,:,:,:,3) = stimulusSpatialModulationOrtho(:,:,:,:,1);
    
    template = stimStruct.template;
    templateOrtho = stimStruct.templateOrtho;
end



function realizedScene = generateScene(theStimulusSpatialModulation,presentationDisplay,stimulusSizeDegs, meanLuminanceCdPerM2)
    realizedScene = sceneFromFile(theStimulusSpatialModulation, 'rgb', ...
                meanLuminanceCdPerM2, presentationDisplay);
    realizedScene = sceneSet(realizedScene, 'fov', stimulusSizeDegs);
end

function d = generateDisplay(viewingDistance, maxLuminance)
    d = displayCreate('LCD-Apple', 'viewing distance', viewingDistance);
    d = displaySet(d, 'gTable', 'linear');
    
    scaleFactor = maxLuminance/displayGet(d,'peak luminance');
    if (scaleFactor > 1)
        spd = displayGet(d,'spd');
        spd = spd * scaleFactor;
        d = displaySet(d, 'spd', spd);
    end
    
    peakLum = displayGet(d,'peak luminance');
    
end

function visualizeSpatialModulations(lowFrequencySpatialModulations, lowFrequencySpatialModulationsOrtho, ...
                    highFrequencySpatialModulations, highFrequencySpatialModulationsOrtho)

    nContrasts = size(lowFrequencySpatialModulations,2);
    
    showAllStimuli = ~true;
    if (showAllStimuli)
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'colsNum', 4, ...
           'rowsNum', nContrasts, ...
           'heightMargin',  0.03, ...
           'widthMargin',    0.03, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.03, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.04);
   
        for theNoiseInstance = 1:noiseInstances
            hFig = figure(1000 + theNoiseInstance); clf;
            set(hFig, 'Position', [10 10 1000 760]);

            for iContrast = 1:nContrasts
            subplot('Position', subplotPosVectors(iContrast,1).v);
            imshow(squeeze(lowFrequencySpatialModulations(theNoiseInstance, iContrast,:,:,:)), [0 1]);
            axis 'image'

            subplot('Position', subplotPosVectors(iContrast,2).v);
            imshow(squeeze(lowFrequencySpatialModulationsOrtho(theNoiseInstance, iContrast,:,:,:)), [0 1]);
            axis 'image'

            subplot('Position', subplotPosVectors(iContrast,3).v);
            imshow(squeeze(highFrequencySpatialModulations(theNoiseInstance, iContrast,:,:,:)), [0 1]);
            axis 'image'

            subplot('Position', subplotPosVectors(iContrast,4).v);
            imshow(squeeze(highFrequencySpatialModulationsOrtho(theNoiseInstance, iContrast,:,:,:)), [0 1]);
            axis 'image'
            end
        end

        colormap(gray(1024));
    end
end