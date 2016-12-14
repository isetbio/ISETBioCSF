function [validationData, extraData] = t_coneCurrentEyeMovementsResponseInstances(varargin)
% [validationData, extraData] = t_coneCurrentEyeMovementsResponseInstances(varargin)
%
% Show how to generate a number of response instances for a given stimulus
% condition.  The default parameters are set up to generate just a single frame
% of the response, but the same tutorial can do temporal sequences with other
% parameter choices.
%
% This tutorial relies on routine
%   colorDetectResponseInstanceFastArrayConstruct
% which does most of the hard work.  The basic principles underlying colorDetectResponseInstanceFastArrayConstruct
% itself is demonstrated in tutorial
%   t_coneCurrentEyeMovementsMovie
% but the actual routine has some tricks to make it go fast.  There is also
% are routine
%   colorDetectResponseInstanceArrayConstruct
% that works more like the tutorial but is slower.
%
% This tutorial saves its output in a .mat file, which cah then read in by
%   t_colorDetectFindPerformance
% which shows how to use the data to find the thresholds.
%
% The returned validation structure allows this routine to be called from a
% validation script driven by the UnitTest toolbox.
%
% The tutorial produces output according to a scheme controlled by the
% specified IBIOColorDetect rwObject.
%
% See also:
%   t_colorGaborRespnseGenerationParams
%   t_coneCurrentEyeMovementsMovie
%	t_colorDetectFindPerformance
%   colorDetectResponseInstanceArrayConstruct
%   colorDetectResponseInstanceFastArrayConstruct
%
% Key/value pairs
%   'rParams' - Value the is the rParams structure to use
%   'testDirectionParams - Value is the testDirectionParams structure to use
%   'emPathType' - Value (one of: 'none', 'frozen', 'random') determines whether we have:
%                  zero eye movements across all trials, 
%                  an emPath that is frozen across all trials,
%                  or a dynamic emPath that changes across trials.
%   'freezeNoise' - true/false (default true).  Freezes all noise so that results are reproducible
%   'compute' - true/false (default true).  Do the computations.
%   'generatePlots' - true/false (default false).  Produce response
%        visualizations.  Set to false when running big jobs on clusters or
%        in parfor loops, as plotting doesn't seem to play well with those
%        conditions.
%   'visualizedResponseNormalization' - How to normalize visualized responses
%        Available options: 'submosaicBasedZscore', 'LMSabsoluteResponseBased', 'LMabsoluteResponseBased', 'MabsoluteResponseBased'
%   'exportPDF' - true/false (default true).  If visualizing responses,
%        export the PDF files.
%   'renderVideo' - true/false (default true).  If visualizing responses, generate
%        the videos.
%   'delete' - true/false (default true).  Delete the response instance
%        files.  Useful for cleaning up big output when we are done with
%        it.  If this is true, output files are deleted at the end.

%% Parse input
p = inputParser;
p.addParameter('rParams',[],@isemptyorstruct);
p.addParameter('testDirectionParams',[],@isemptyorstruct);
p.addParameter('emPathType','none',@ischar);
p.addParameter('freezeNoise',true,@islogical);
p.addParameter('compute',true,@islogical);
p.addParameter('generatePlots',false,@islogical);
p.addParameter('exportPDF',true,@islogical);
p.addParameter('renderVideo',true,@islogical);
p.addParameter('delete',false',@islogical);
p.addParameter('visualizedResponseNormalization', 'submosaicBasedZscore', @ischar);
p.parse(varargin{:});
rParams = p.Results.rParams;
testDirectionParams = p.Results.testDirectionParams;

%% Clear
if (nargin == 0)
    ieInit; close all;
end

%% Get the parameters we need
%
% t_colorGaborResponseGenerationParams returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
if (isempty(rParams))
    rParams = responseParamsGenerate;
    
    % Override some defult parameters
    %
    % Set duration equal to sampling interval to do just one frame.
    rParams.temporalParams.stimulusDurationInSeconds = 200/1000;
    rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.stimulusDurationInSeconds;
    rParams.temporalParams.secondsToInclude = rParams.temporalParams.stimulusDurationInSeconds;
    
    rParams.mosaicParams.integrationTimeInSeconds = rParams.temporalParams.stimulusDurationInSeconds;
    rParams.mosaicParams.isomerizationNoise = 'random';         % Type coneMosaic.validNoiseFlags to get valid values
    rParams.mosaicParams.osNoise = 'random';                    % Type outerSegment.validNoiseFlags to get valid values
    rParams.mosaicParams.osModel = 'Linear';
end


% Set the emPathType
if (~isempty(p.Results.emPathType))
    rParams.temporalParams.emPathType = p.Results.emPathType;
end
    
% Fix random number generator so we can validate output exactly
if (p.Results.freezeNoise)
     fprintf(2, '\n%s: freezing all noise \n', mfilename);
     rng(1);
     if (strcmp(rParams.mosaicParams.isomerizationNoise, 'random'))
         fprintf(2, '\tmosaicParams.isomerizationNoise was set to ''%s'', setting it to ''frozen''.\n', rParams.mosaicParams.isomerizationNoise);
         rParams.mosaicParams.isomerizationNoise = 'frozen';
     end
     if (strcmp(rParams.mosaicParams.osNoise, 'random'))
         fprintf(2, '\tmosaicParams.osNoise was set to ''%s'', setting it to ''frozen''.\n', rParams.mosaicParams.osNoise);
         rParams.mosaicParams.osNoise = 'frozen';
     end
end

%% Parameters that define the LM instances we'll generate here
if (isempty(testDirectionParams))
    testDirectionParams = instanceParamsGenerate;
end

% The constant params list
constantParamsList = {rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams,  rParams.backgroundParams, testDirectionParams};


%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
theProgram = mfilename;

%% The computing happens here, if we are doing it
if (p.Results.compute)
    % Create the optics
    theOI = colorDetectOpticalImageConstruct(rParams.oiParams);
    
    % Create the cone mosaic
    theMosaic = colorDetectConeMosaicConstruct(rParams.mosaicParams);
    
    %% Define color modulation list
    switch (testDirectionParams.instanceType)
        case 'LMPlane'
            % Directions
            testConeContrasts = testConeContrastsFromTestDirectionParams(rParams,testDirectionParams);
            
            % Contrasts
            if (strcmp(testDirectionParams.contrastScale, 'linear'))
                testContrasts = linspace(testDirectionParams.lowContrast, testDirectionParams.highContrast, testDirectionParams.nContrastsPerDirection);
            else
                testContrasts = logspace(log10(testDirectionParams.lowContrast), log10(testDirectionParams.highContrast), testDirectionParams.nContrastsPerDirection);
            end
            
        case 'LMSPlane'
            % Directions
            testConeContrasts = testConeContrastsFromTestDirectionParams(rParams,testDirectionParams);
            
            % Contrasts
            if (strcmp(testDirectionParams.contrastScale, 'linear'))
                testContrasts = linspace(testDirectionParams.lowContrast, testDirectionParams.highContrast, testDirectionParams.nContrastsPerDirection);
            else
                testContrasts = logspace(log10(testDirectionParams.lowContrast), log10(testDirectionParams.highContrast), testDirectionParams.nContrastsPerDirection);
            end
            
            
        case 'contrasts'
            % Contrasts
            if (strcmp(testDirectionParams.contrastScale, 'linear'))
                testContrasts = linspace(testDirectionParams.lowContrast, testDirectionParams.highContrast, testDirectionParams.nContrastsPerDirection);
            else
                testContrasts = logspace(log10(testDirectionParams.lowContrast), log10(testDirectionParams.highContrast), testDirectionParams.nContrastsPerDirection);
            end
            
            % Set up parfor condition structs. Just dummy up a cone
            % contrasts argument for this call, we won't use the returned
            % contrasts for this case below.
            testConeContrasts = NaN*zeros(3,1);
        otherwise
            error('Unknown test instance type passed');
    end
    parforConditionStructs = responseGenerationParforConditionStructsGenerate(testConeContrasts,testContrasts);
    nParforConditions = length(parforConditionStructs);
    parforRanSeeds = randi(1000000,nParforConditions,1)+1;

    % Generate data for the no stimulus condition
    colorModulationParamsNull = rParams.colorModulationParams;
    colorModulationParamsNull.coneContrasts = [0 0 0]';
    colorModulationParamsNull.contrast = 0;
    stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.2f', ...
        colorModulationParamsNull.coneContrasts(1), colorModulationParamsNull.coneContrasts(2), colorModulationParamsNull.coneContrasts(3), colorModulationParamsNull.contrast);
    [responseInstanceArray,noiseFreeIsomerizations, noiseFreePhotocurrents] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, testDirectionParams.trialsNum, ...
        rParams.spatialParams, rParams.backgroundParams, colorModulationParamsNull, rParams.temporalParams, theOI, theMosaic, 'seed', 1);
 
    noStimData = struct(...
        'testContrast', colorModulationParamsNull.contrast, ...
        'testConeContrasts', colorModulationParamsNull.coneContrasts, ...
        'stimulusLabel', stimulusLabel, ...
        'responseInstanceArray',responseInstanceArray, ...
        'noiseFreeIsomerizations',noiseFreeIsomerizations, ...
        'noiseFreePhotocurrents', noiseFreePhotocurrents);
    
    % Write the no cone contrast data and some extra facts we need
    paramsList = constantParamsList;
    paramsList{numel(paramsList)+1} = colorModulationParamsNull;
    rwObject.write('responseInstances',noStimData,paramsList,theProgram);
      
    % Save the other data we need for use by the classifier preprocessing subroutine
    ancillaryData = struct(...
        'testConeContrasts', testConeContrasts, ...
        'testContrasts', testContrasts, ...
        'rParams', rParams, ...
        'instanceParams', testDirectionParams);
    ancillaryData.parforConditionStructs = parforConditionStructs;
    rwObject.write('ancillaryData',ancillaryData,paramsList,theProgram);
    
    %% Generate data for all the examined stimuli
    %
    % It is possible that the parfor loop will not work for you, depending
    % on your Matlab configuration.  In this case, change it to a for loop.
   
    % Loop over color directions
    %
    % Note tha as the mosaic handle object (and any other handle objects)
    % enter the parfor loop, a copy local to each worker is created.  This
    % is the behavior we want, so that the isomerizations calculations for
    % each loop iteration don't step on each other.  Also note that any
    % changes to the mosaic object inside the loop do not get propagted
    % back out -- the mosaic object itself is the same at the end of the
    % parfor as at the start.  This is also OK here, but might be confusing
    % under other circumstances.
    tic;
    stimDataForValidation = cell(nParforConditions,1);

    parfor kk = 1:nParforConditions  
        thisConditionStruct = parforConditionStructs{kk};
        colorModulationParamsTemp = rParams.colorModulationParams;
        colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
        colorModulationParamsTemp.contrast = thisConditionStruct.contrast;
        
        % Make noisy instances for each contrast
        stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.5f',...
            colorModulationParamsTemp.coneContrasts(1), colorModulationParamsTemp.coneContrasts(2), colorModulationParamsTemp.coneContrasts(3), colorModulationParamsTemp.contrast);
        [responseInstanceArray,noiseFreeIsomerizations, noiseFreePhotocurrents] = ...
            colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, testDirectionParams.trialsNum, ...
                rParams.spatialParams, rParams.backgroundParams, colorModulationParamsTemp, rParams.temporalParams, theOI, theMosaic, 'seed', parforRanSeeds(kk));
        stimData = struct(...
            'testContrast', colorModulationParamsTemp.contrast, ...
            'testConeContrasts', colorModulationParamsTemp.coneContrasts, ...
            'stimulusLabel', stimulusLabel, ...
            'responseInstanceArray',responseInstanceArray, ...
            'noiseFreeIsomerizations',noiseFreeIsomerizations, ...
            'noiseFreePhotocurrents', noiseFreePhotocurrents);
        
        % Save some data for validation in first loop
        if (kk == 1)
            s = struct();
            savedTrial = 1;
            s.theMosaicIsomerizations = squeeze(responseInstanceArray.theMosaicIsomerizations(savedTrial,:,:,:));
            if (isempty(responseInstanceArray.theMosaicPhotocurrents))
                s.theMosaicPhotocurrents = [];
            else
                s.theMosaicPhotocurrents = squeeze(responseInstanceArray.theMosaicPhotocurrents(savedTrial,:,:,:));
            end
            s.theMosaicEyeMovements = squeeze(responseInstanceArray.theMosaicEyeMovements(savedTrial,:,:));
            s.timeAxis = responseInstanceArray.timeAxis;
            s.photocurrentTimeAxis = responseInstanceArray.photocurrentTimeAxis;
            stimDataForValidation{kk} = s;
        end
        
        % Save data for this color direction/contrast pair
        paramsList = constantParamsList;
        paramsList{numel(paramsList)+1} = colorModulationParamsTemp;
        rwObject.write('responseInstances',stimData,paramsList,theProgram);
    end
    fprintf('Finished generating responses in %2.2f minutes\n', toc/60);
    
    %% Validation data
    %
    % Since we want to store some stimData and that isn't saved outside the
    % parfor loop (and is in fact a bit hard to save outside the parfor
    % loop)
    if (nargout > 0)
        savedTrial = 1;
        noStimValidationData.theMosaicIsomerizations = squeeze(noStimData.responseInstanceArray.theMosaicIsomerizations(savedTrial,:,:,:));
        if (isempty(noStimData.responseInstanceArray.theMosaicPhotocurrents))
            noStimValidationData.theMosaicPhotocurrents = [];
        else
            noStimValidationData.theMosaicPhotocurrents = squeeze(noStimData.responseInstanceArray.theMosaicPhotocurrents(savedTrial,:,:,:));
        end
        noStimValidationData.theMosaicEyeMovements = squeeze(noStimData.responseInstanceArray.theMosaicEyeMovements(savedTrial,:,:));
        noStimValidationData.timeAxis = noStimData.responseInstanceArray.timeAxis;
        noStimValidationData.photocurrentTimeAxis = noStimData.responseInstanceArray.photocurrentTimeAxis;
        
        validationData.noStimData = noStimValidationData;
        validationData.stimData = stimDataForValidation;
        extraData.ancillaryData = ancillaryData;
        extraData.p.Results = p.Results;
    end
end

%% Visualize
if (p.Results.generatePlots)

    % How many istances to visualize
    instancesToVisualize = 1:5;
    
    % Load data
    paramsList = constantParamsList;
    paramsList{numel(paramsList)+1} = colorModulationParamsNull;    
    noStimData = rwObject.read('responseInstances',paramsList,theProgram);
    ancillaryData = rwObject.read('ancillaryData',paramsList,theProgram);
    
    rParams = ancillaryData.rParams;
    parforConditionStructs = ancillaryData.parforConditionStructs;
    nParforConditions = length(parforConditionStructs); 
    for kk = 1:nParforConditions 
         thisConditionStruct = parforConditionStructs{kk};
         colorModulationParamsTemp = rParams.colorModulationParams;
         colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
         colorModulationParamsTemp.contrast = thisConditionStruct.contrast;

         paramsList = constantParamsList;
         paramsList{numel(paramsList)+1} = colorModulationParamsTemp;    
         stimData = rwObject.read('responseInstances',paramsList,theProgram);
         visualizeResponses(theMosaic, stimData, noStimData, p.Results.visualizedResponseNormalization, kk, nParforConditions, instancesToVisualize);
    end
end


%% Delete response instance files.
%
% Doesn't delete figures or parent directories.
if (p.Results.delete)
    rwObject.delete('responseInstances',paramsList,theProgram);
    rwObject.delete('ancillaryData',paramsList,theProgram);
    for kk = 1:nParforConditions 
        thisConditionStruct = parforConditionStructs{kk};
        colorModulationParamsTemp = rParams.colorModulationParams;
        colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
        colorModulationParamsTemp.contrast = thisConditionStruct.contrast;
        paramsList = constantParamsList;
        paramsList{numel(paramsList)+1} = colorModulationParamsTemp;
    
       % paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
        rwObject.delete('responseInstances',paramsList,theProgram);   
    end   
end
end


function visualizeResponses(theMosaic, stimData, noStimData, responseNormalization, condIndex, condsNum, instancesToVisualize)
         
    instancesNum = size(stimData.responseInstanceArray.theMosaicIsomerizations,1);
    if (instancesNum < 1)
        return;
    end
    
    for coneIndex = 2:4
        submosaicConeIndices{coneIndex} = find(theMosaic.pattern==coneIndex);
    end

    timeBins = numel(stimData.responseInstanceArray.timeAxis);
    if (timeBins == 1)
        if (isa(theMosaic, 'coneMosaicHex'))
            coneDims = 2;
        else
            coneDims = [2 3];
        end
    elseif (timeBins > 1)
        if (isa(theMosaic, 'coneMosaicHex'))
            coneDims = 2;
        else
            coneDims = [2 3];
        end
    else
        error('timeBins = %d', timeBins)
    end
    
    photocurrents = [];
    if (strcmp(responseNormalization, 'LMSabsoluteResponseBased')) || (strcmp(responseNormalization, 'LMabsoluteResponseBased')) || (strcmp(responseNormalization, 'MabsoluteResponseBased')) 
        % Max from L- and M-cone mosaics
        if (strcmp(responseNormalization, 'LMSabsoluteResponseBased')) 
            normalizationConeIndices = [submosaicConeIndices{2}; submosaicConeIndices{3}; submosaicConeIndices{4}];
        elseif (strcmp(responseNormalization, 'LMabsoluteResponseBased')) 
            normalizationConeIndices = [submosaicConeIndices{2}; submosaicConeIndices{3};];
        elseif (strcmp(responseNormalization, 'MabsoluteResponseBased')) 
            % Max from M-cone mosaic only
            normalizationConeIndices = [submosaicConeIndices{3}];
        else
            error('unknown normalization: ''%s''.', responseNormalization);
        end
        [absorptions, minAbsorptions, maxAbsorptions] = coneIndicesBasedScaling(stimData.responseInstanceArray.theMosaicIsomerizations, coneDims, normalizationConeIndices, true);
        [noiseFreeIsomerizations, minNoiseFreeIsomerizations, maxNoiseFreeIsomerizations] = coneIndicesBasedScaling(stimData.noiseFreeIsomerizations, coneDims, normalizationConeIndices, false);
        if (~isempty(stimData.noiseFreePhotocurrents))
            [photocurrents, minPhotocurrents, maxPhotocurrents] = coneIndicesBasedScaling(stimData.responseInstanceArray.theMosaicPhotocurrents, coneDims, normalizationConeIndices, true);
            [noiseFreePhotocurrents, minNoiseFreePhotocurrents, maxNoiseFreePhotocurrents] = coneIndicesBasedScaling(stimData.noiseFreePhotocurrents, coneDims, normalizationConeIndices, false);
        end
    elseif strcmp(responseNormalization, 'submosaicBasedZscore')
         [absorptions, minAbsorptions, maxAbsorptions] = submosaicBasedZscore(stimData.responseInstanceArray.theMosaicIsomerizations, noStimData.responseInstanceArray.theMosaicIsomerizations, coneDims, submosaicConeIndices, true);
         [noiseFreeIsomerizations, minNoiseFreeIsomerizations, maxNoiseFreeIsomerizations] = submosaicBasedZscore(stimData.noiseFreeIsomerizations, noStimData.noiseFreeIsomerizations,coneDims, submosaicConeIndices, false);
         if (~isempty(stimData.noiseFreePhotocurrents))
             [photocurrents, minPhotocurrents, maxPhotocurrents] = submosaicBasedZscore(stimData.responseInstanceArray.theMosaicPhotocurrents, noStimData.responseInstanceArray.theMosaicPhotocurrents, coneDims,  submosaicConeIndices, true);
             [noiseFreePhotocurrents, minNoiseFreePhotocurrents, maxNoiseFreePhotocurrents] = submosaicBasedZscore(stimData.noiseFreePhotocurrents, noStimData.noiseFreePhotocurrents, coneDims,  submosaicConeIndices, false);
         end
    else
        error('Unknown responseNormalization method: ''%s''.', responseNormalization);
    end

     
    if (isa(theMosaic, 'coneMosaicHex'))
        fprintf(2, '\nConeHex visualization not implemented yet\n');
        return;
    end
    
    absorptionsTimeAxis = stimData.responseInstanceArray.timeAxis;
    photocurrentsTimeAxis = stimData.responseInstanceArray.photocurrentTimeAxis;
         
    hFig = figure(100+condIndex); clf;
    set(hFig, 'Position', [10 10 1100 1050], 'Color', [1 1 1], 'Name', stimData.stimulusLabel);

    if (isempty(photocurrents))
        subplotRows = 1;
    else
        subplotRows = 2;
    end
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', subplotRows, ...
           'colsNum', 2, ...
           'heightMargin',   0.06, ...
           'widthMargin',    0.13, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.08, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.03);
       
    
    for instanceIndex = instancesToVisualize
         for tBin = 1: numel(absorptionsTimeAxis)  
             
            % Instance absorptions on the left
            subplot('Position', subplotPosVectors(1,1).v);
            imagesc(squeeze(absorptions(instanceIndex, :,:,tBin))); axis 'image'
            set(gca, 'CLim', [0 1], 'FontSize', 14);
            title(sprintf('absorptions\ninstance %d/%d (t: %2.1fms)', instanceIndex, instancesNum, absorptionsTimeAxis(tBin)*1000));

            % Add colorbar
            originalPosition = get(gca, 'position');

            ticks = 0:0.2:1.0;
            delta = (maxAbsorptions-minAbsorptions)*0.2;
            if strcmp(responseNormalization, 'submosaicBasedZscore')
                tickLabels = (minAbsorptions:delta:maxAbsorptions);
                colorbarLabel = sprintf('absorptions z-score (%2.2fms)', theMosaic.integrationTime*1000);
            else
                tickLabels = (minAbsorptions:delta:maxAbsorptions);
                colorbarLabel = sprintf('absorptions (R*/cone/%2.2fms)', theMosaic.integrationTime*1000);
            end
            hCbar = colorbar('Ticks', ticks, 'TickLabels', sprintf('%2.2f\n',tickLabels));
            hCbar.Orientation = 'vertical'; 
            hCbar.Label.String = colorbarLabel;
            hCbar.FontSize = 14; 
            hCbar.FontName = 'Menlo'; 
            hCbar.FontWeight = 'Bold'; 
            hCbar.Color = [0.2 0.2 0.2];
            % The addition changes the figure size, so undo this change
            newPosition = get(gca, 'position');
            set(gca,'position',[newPosition(1) newPosition(2) originalPosition(3) originalPosition(4)]);
            
            
            % Noise-free isomerizations on the right
            subplot('Position', subplotPosVectors(1,2).v);
            imagesc(squeeze(noiseFreeIsomerizations(:,:,tBin))); axis 'image'
            set(gca, 'CLim', [0 1], 'FontSize', 14);
            title(sprintf('noise-free absorptions\n cond: %d/%d  (t: %2.1fms)', condIndex, condsNum,  absorptionsTimeAxis(tBin)*1000));

            % Add colorbar
            originalPosition = get(gca, 'position');

            ticks = 0:0.2:1.0;
            delta = (maxNoiseFreeIsomerizations-minNoiseFreeIsomerizations)*0.2;
            tickLabels = (minNoiseFreeIsomerizations:delta:maxNoiseFreeIsomerizations);
            if strcmp(responseNormalization, 'submosaicBasedZscore')
                colorbarLabel = sprintf('noise-free absorptions zscore (%2.2fms)', theMosaic.integrationTime*1000);
            else
                colorbarLabel = sprintf('noise-free absorptions (R*/cone/%2.2fms)', theMosaic.integrationTime*1000);
            end
            hCbar = colorbar('Ticks', ticks, 'TickLabels', sprintf('%2.2f\n',tickLabels));
            hCbar.Orientation = 'vertical'; 
            hCbar.Label.String = colorbarLabel;
            hCbar.FontSize = 14; 
            hCbar.FontName = 'Menlo'; 
            hCbar.FontWeight = 'Bold'; 
            hCbar.Color = [0.2 0.2 0.2];
            % The addition changes the figure size, so undo this change
            newPosition = get(gca, 'position');
            set(gca,'position',[newPosition(1) newPosition(2) originalPosition(3) originalPosition(4)]);
            

            if (~isempty(photocurrents))
                % Instance photocurrents on the left
                subplot('Position', subplotPosVectors(2,1).v);
                imagesc(squeeze(photocurrents(instanceIndex,:,:,tBin))); axis 'image'
                set(gca, 'CLim', [0 1], 'FontSize', 14);
                title(sprintf('Photocurrents\n instance %d/%d (t: %2.1fms)', instanceIndex, instancesNum, photocurrentsTimeAxis(tBin)*1000));

                % Add colorbar
                originalPosition = get(gca, 'position');
                ticks = 0:0.2:1.0;
                delta = (maxPhotocurrents-minPhotocurrents)*0.2;
                if strcmp(responseNormalization, 'submosaicBasedZscore')
                    tickLabels = minPhotocurrents:delta:maxPhotocurrents;
                    colorbarLabel = sprintf('photocurrents z-score');
                else
                    tickLabels = minPhotocurrents:delta:maxPhotocurrents;
                    colorbarLabel = sprintf('photocurrents (pAmps)');
                end
                hCbar = colorbar('Ticks', ticks, 'TickLabels', sprintf('%2.1f\n',tickLabels));
                hCbar.Orientation = 'vertical'; 
                hCbar.Label.String = colorbarLabel;
                hCbar.FontSize = 14; 
                hCbar.FontName = 'Menlo'; 
                hCbar.FontWeight = 'Bold'; 
                hCbar.Color = [0.2 0.2 0.2];
                % The addition changes the figure size, so undo this change
                newPosition = get(gca, 'position');
                set(gca,'position',[newPosition(1) newPosition(2) originalPosition(3) originalPosition(4)]);
                
                
                % Mean photocurrents on the right
                subplot('Position', subplotPosVectors(2,2).v);
                imagesc(squeeze(noiseFreePhotocurrents(:,:,tBin))); axis 'image'
                set(gca, 'CLim', [0 1], 'FontSize', 14);
                title(sprintf('noise-free photocurrents \n cond: %d/%d  (t: %2.1fms)', condIndex, condsNum, photocurrentsTimeAxis(tBin)*1000));
                
                % Add colorbar
                originalPosition = get(gca, 'position');
                ticks = 0:0.2:1.0;
                delta = (maxNoiseFreePhotocurrents-minNoiseFreePhotocurrents)*0.2;
                if strcmp(responseNormalization, 'submosaicBasedZscore')
                    tickLabels = minNoiseFreePhotocurrents:delta:maxNoiseFreePhotocurrents;
                    colorbarLabel = sprintf('photocurrents z-score');
                else
                    tickLabels = minNoiseFreePhotocurrents:delta:maxNoiseFreePhotocurrents;
                    colorbarLabel = sprintf('photocurrents (pAmps)');
                end
                hCbar = colorbar('Ticks', ticks, 'TickLabels', sprintf('%2.1f\n',tickLabels));
                hCbar.Orientation = 'vertical'; 
                hCbar.Label.String = colorbarLabel;
                hCbar.FontSize = 14; 
                hCbar.FontName = 'Menlo'; 
                hCbar.FontWeight = 'Bold'; 
                hCbar.Color = [0.2 0.2 0.2];
                % The addition changes the figure size, so undo this change
                newPosition = get(gca, 'position');
                set(gca,'position',[newPosition(1) newPosition(2) originalPosition(3) originalPosition(4)]);
            end
            
            colormap(bone(1024));
            drawnow;
         end
     end
end 
    

function [responseDataZscore, minZscore, maxZscore] = submosaicBasedZscore(responseData, noResponseData, coneDims, submosaicConeIndices, isInstanceData)
    if (numel(coneDims) == 2)
        originalResponseDataDims = size(responseData);
        if (isInstanceData)
            responseData = reshape(responseData, [size(responseData,1) size(responseData,2)*size(responseData,3) size(responseData,4)]);
            noResponseData = reshape(noResponseData, [size(noResponseData,1) size(noResponseData,2)*size(noResponseData,3) size(noResponseData,4)]);
        else
            responseData = reshape(responseData, [size(responseData,1)*size(responseData,2) size(responseData,3)]);
            noResponseData = reshape(noResponseData, [size(noResponseData,1)*size(noResponseData,2) size(noResponseData,3)]);
        end
    end
    
    responseDataZscore = responseData;
    for coneIndex = 2:4
        if (isInstanceData)
            subMosaicResponseData = responseData(:, submosaicConeIndices{coneIndex},:);
            subMosaicNoResponseData = noResponseData(:, submosaicConeIndices{coneIndex},:);
        else
            % noise-free data
            subMosaicResponseData = responseData(submosaicConeIndices{coneIndex},:);
            subMosaicNoResponseData = noResponseData(submosaicConeIndices{coneIndex},:);
        end
        
        % subtract mean over all cones of a particular type (across all instances and time) for the noResponse data
        meanSubMosaicNoResponse = mean(subMosaicNoResponseData(:));
        subMosaicResponseData = (subMosaicResponseData - meanSubMosaicNoResponse);
        
        % divide by std over all cones of a particular type (across all instances and time) for the noResponse data
        if (isInstanceData)
            stdSubMosaicNoResponse = std(subMosaicNoResponseData(:)); 
            subMosaicResponseData = subMosaicResponseData/stdSubMosaicNoResponse;
        end
        
        if (isInstanceData)
            responseDataZscore(:, submosaicConeIndices{coneIndex},:) = subMosaicResponseData;
        else
            % noise-free data
            responseDataZscore(submosaicConeIndices{coneIndex},:) = subMosaicResponseData;
        end
    end % coneIndex

    % Normalize
    if (isInstanceData)
        maxZscore = 1.0;
        minZscore = -maxZscore;
    else
        maxZscore = max(abs(responseDataZscore(:)));
        minZscore = -maxZscore;
    end
    responseDataZscore = (responseDataZscore-minZscore)/(maxZscore-minZscore);  
    
    % Back to original shape
    if (numel(coneDims) == 2)
        responseDataZscore = reshape(responseDataZscore, originalResponseDataDims);
    end
end
 
 
 function [scaledResp, minResp, maxResp] = coneIndicesBasedScaling(responseData, coneDims, coneIndices, isInstanceData)
    if (numel(coneDims) == 2)
        originalResponseDataDims = size(responseData);
        if (isInstanceData)
            responseData = reshape(responseData, [size(responseData,1) size(responseData,2)*size(responseData,3) size(responseData,4)]);
        else
            responseData = reshape(responseData, [size(responseData,1)*size(responseData,2) size(responseData,3)]);
        end
    end
    
    if (isInstanceData)
        subMosaicResponseData = responseData(:, coneIndices,:);
    else
        % noise-free data
        subMosaicResponseData = responseData(coneIndices,:);
    end
        
    maxResp = max(subMosaicResponseData(:));
    minResp = min(subMosaicResponseData(:));
    scaledResp = (responseData-minResp)/(maxResp-minResp);
    
    % Back to original shape
    if (numel(coneDims) == 2)
        scaledResp = reshape(scaledResp, originalResponseDataDims);
    end
 end
 
