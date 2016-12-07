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
%   'emPathType' - Value (one of: 'Zero', 'Frozen', 'Dynamic') determines whether we have:
%                  zero eye movements across all trials, 
%                  an emPath that is frozen across all trials,
%                  or a dynamic emPath that changes across trials.
%   'freezeNoise' - true/false (default true).  Freezes all noise so that results are reproducible
%   'compute' - true/false (default true).  Do the computations.
%   'generatePlots' - true/false (default false).  Produce response
%        visualizations.  Set to false when running big jobs on clusters or
%        in parfor loops, as plotting doesn't seem to play well with those
%        conditions.
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
    rParams.temporalParams.simulationTimeStepSecs = 200/1000;
    rParams.temporalParams.stimulusDurationInSeconds = rParams.temporalParams.simulationTimeStepSecs;
    rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.simulationTimeStepSecs;
    rParams.temporalParams.secondsToInclude = rParams.temporalParams.simulationTimeStepSecs;
    
    rParams.mosaicParams.timeStepInSeconds = rParams.temporalParams.simulationTimeStepSecs;
    rParams.mosaicParams.integrationTimeInSeconds = rParams.mosaicParams.timeStepInSeconds;
    rParams.mosaicParams.isomerizationNoise = 'random';         % select from {'random', 'frozen', 'none'}
    rParams.mosaicParams.osNoise = 'random';                    % select from {'random', 'frozen', 'none'}
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

%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
theProgram = mfilename;

%% The computing happens here, if we are doing it
if (p.Results.compute)
    % Create the optics
    theOI = colorDetectOpticalImageConstruct(rParams.oiParams);
    
    % Create the cone mosaic
    rParams.mosaicParams.fieldOfViewDegs = rParams.spatialParams.fieldOfViewDegs;
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
    colorModulationParamsTemp = rParams.colorModulationParams;
    colorModulationParamsTemp.coneContrasts = [0 0 0]';
    colorModulationParamsTemp.contrast = 0;
    stimulusLabel = sprintf('LMS=%2.2f,%2.2f,%2.2f,Contrast=%2.2f', ...
        colorModulationParamsTemp.coneContrasts(1), colorModulationParamsTemp.coneContrasts(2), colorModulationParamsTemp.coneContrasts(3), colorModulationParamsTemp.contrast);
    [responseInstanceArray,noiseFreeIsomerizations, noiseFreePhotocurrents] = colorDetectResponseInstanceArrayFastConstruct(stimulusLabel, testDirectionParams.trialsNum, ...
        rParams.spatialParams, rParams.backgroundParams, colorModulationParamsTemp, rParams.temporalParams, theOI, theMosaic, 'seed', 1);
 
    noStimData = struct(...
        'testContrast', colorModulationParamsTemp.contrast, ...
        'testConeContrasts', colorModulationParamsTemp.coneContrasts, ...
        'stimulusLabel', stimulusLabel, ...
        'responseInstanceArray',responseInstanceArray, ...
        'noiseFreeIsomerizations',noiseFreeIsomerizations, ...
        'noiseFreePhotocurrents', noiseFreePhotocurrents);
    
    % Write the no cone contrast data and some extra facts we need
    paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
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
    %rState = rng;
    parfor kk = 1:nParforConditions          
        %rng(parforRanSeeds(kk));
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
            stimDataForValidation{kk} = stimData.responseInstanceArray(1);
        end
        
        % Save data for this color direction/contrast pair
        paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
        rwObject.write('responseInstances',stimData,paramsList,theProgram);
    end
    %rng(rState);
    fprintf('Finished generating responses in %2.2f minutes\n', toc/60);
    
    %% Validation data
    %
    % Since we want to store some stimData and that isn't saved outside the
    % parfor loop (and is in fact a bit hard to save outside the parfor
    % loop)
    if (nargout > 0)
        validationData.noStimData = noStimData.responseInstanceArray(1);
        validationData.stimData = stimDataForValidation{1};
        extraData.ancillaryData = ancillaryData;
        extraData.p.Results = p.Results;
    end
end

%% Visualize
%
% This is not yet implemented, but here is where the code would go.  And, a
% very early version of the visualization is in the commented out code
% immediately below
if (p.Results.generatePlots)

    noStimData = rwObject.read('responseInstances',paramsList,theProgram);
    ancillaryData = rwObject.read('ancillaryData',paramsList,theProgram);
    parforConditionStructs = ancillaryData.parforConditionStructs;

    rParams = ancillaryData.rParams;
    nParforConditions = length(parforConditionStructs); 
    for kk = 1:nParforConditions 
         thisConditionStruct = parforConditionStructs{kk};
         colorModulationParamsTemp = rParams.colorModulationParams;
         colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
         colorModulationParamsTemp.contrast = thisConditionStruct.contrast;
         paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
         stimData = rwObject.read('responseInstances',paramsList,theProgram);   
         instancesNum = size(stimData.responseInstanceArray,1);
         
         hFig = figure(100); clf;
         set(hFig, 'Position', [10 10 1400 800], 'Color', [1 1 1]);
         
         size(stimData.noiseFreeIsomerizations)
         
         % normalize submosaic responses with respect to their mean
         responseNormalization = 'absoluteLevel';
         if strcmp(responseNormalization, 'absoluteLevel')
            maxNoiseFreeIsomerization = max(stimData.noiseFreeIsomerizations(:));
            minNoiseFreeIsomerization = min(stimData.noiseFreeIsomerizations(:));
            margin = 0.1*maxNoiseFreeIsomerization;
            maxNoiseFreeIsomerization = maxNoiseFreeIsomerization + margin;
            minNoiseFreeIsomerization = minNoiseFreeIsomerization - margin;
         end
         
         noiseFreeAbsorptions = (stimData.noiseFreeIsomerizations-minNoiseFreeIsomerization)/(maxNoiseFreeIsomerization-minNoiseFreeIsomerization);
                 
         for instanceIndex = 1:instancesNum     
             data = stimData.responseInstanceArray(instanceIndex);
             absorptions = data.theMosaicIsomerizations;
             photocurrents = data.theMosaicPhotoCurrents;
             absorptionsTimeAxis = data.timeAxis;
             photocurrentsTimeAxis = data.photocurrentTimeAxis;
             
             if (numel(absorptionsTimeAxis) == 1)
                 absorptions = (absorptions-minNoiseFreeIsomerization)/(maxNoiseFreeIsomerization-minNoiseFreeIsomerization);
                 
                 subplot(1,2,1)
                 imagesc(absorptions); axis 'image'
                 set(gca, 'CLim', [0 1], 'FontSize', 14);
                 title(sprintf('instance:%d/%d', instanceIndex, instancesNum));
                 
                 subplot(1,2,2)
                 imagesc(noiseFreeAbsorptions); axis 'image'
                 set(gca, 'CLim', [0 1], 'FontSize', 14);
                 title(sprintf('cond: %d/%d %s', kk, nParforConditions, stimData.stimulusLabel));
                 
                 % Add colorbar
                 originalPosition = get(gca, 'position');
                 
                 ticks = 0:0.2:1.0;
                 delta = (maxNoiseFreeIsomerization-minNoiseFreeIsomerization)*0.2;
                 tickLabels = minNoiseFreeIsomerization:delta:maxNoiseFreeIsomerization;
                 hCbar = colorbar('Ticks', ticks, 'TickLabels', sprintf('%2.1f\n',tickLabels));
                 hCbar.Orientation = 'vertical'; 
                 hCbar.Label.String = sprintf('absorptions (R*/cone/%2.2fms)', theMosaic.integrationTime*1000);
                 hCbar.FontSize = 14; 
                 hCbar.FontName = 'Menlo'; 
                 hCbar.Color = [0.2 0.2 0.2];
                 % The addition changes the figure size, so undo this change
                 newPosition = get(gca, 'position');
                 set(gca,'position',[newPosition(1) newPosition(2) originalPosition(3) originalPosition(4)]);
        
                 colormap(gray);
                 drawnow;
             end
         end
     end   
end

% % THIS IS BROKEN AND NEES TO BE UPDATED TO NEW DATA FORMAT AND rwObject
% % LAND.
% %
% % Also, the time numbers on the videos do not seem to correspond to the
% % stimulus peak at time 0, and the videos look screwy in the no noise
% % condition.
% if (p.Results.generatePlots)
    % fprintf('\nVisualizing responses ...\n');
    % for ii = 1:size(testConeContrasts,2)
    %     for jj = 1:numel(testContrasts)
    %         stimulusLabel = sprintf('LMS_%2.2f_%2.2f_%2.2f_Contrast_%2.2f', testConeContrasts(1,ii), testConeContrasts(2,ii), testConeContrasts(3,ii), testContrasts(jj));
    %         s = theStimData{ii, jj};
    % 
    %         % Visualize a few response instances only
    %         for iTrial = 1:2
    %             figHandle = visualizeResponseInstance(conditionDir, s.responseInstanceArray(iTrial), stimulusLabel, theMosaic, iTrial, testDirectionParams.trialsNum, p.Results.renderVideo);
    %             if (p.Results.exportPDF)
    %                 figFileNames{ii, jj, iTrial} = ...
    %                     fullfile(colorGaborDetectOutputDir(conditionDir),sprintf('%s_Trial%dOf%d.pdf', stimulusLabel, iTrial, testDirectionParams.trialsNum),'figures');
    %                 NicePlot.exportFigToPDF(figFileNames{ii, jj, iTrial}, figHandle, 300);
    %             end
    %         end
    %     end
    % end
    % 
    % % Export summary PDF with all responses
    % if (p.Results.exportPDF)
    %     summaryPDF = fullfile(colorGaborDetectOutputDir(conditionDir,'figures'), 'AllInstances.pdf');
    %     fprintf('Exporting a summary PDF with all response instances in %s\n', summaryPDF);
    %     NicePlot.combinePDFfilesInSinglePDF(figFileNames(:), summaryPDF);
    % end
% end

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
        paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
        rwObject.delete('responseInstances',paramsList,theProgram);   
    end   
end
end


