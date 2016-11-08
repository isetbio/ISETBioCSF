function c_PoirsonAndWandell96ReplicateParfor
% c_PoirsonAndWandell96Replicate
%
% Compute color detection thresholds to replicate the Poirson & Wandell 1996
         
    close all
    
    % Whether to display all the params
    paramsVerbosity = 0;
    
    % Whether to save videos and stills
    exportLMSresponseTraceStills = true;
    exportMosaic2DActivationStillsAndVideo = true;
    
    % Set to true if cone mosaic does not exist
    recomputeConeMosaic = false;
    
    % How many response instances to generate
    instancesNum = 1000; 
    
    % How much data to save for classification
    temporalResponseToIncludeInUnitsOfIntegrationTime = 0;     % Only peak response
    temporalResponseToIncludeInUnitsOfIntegrationTime = 17.0;  
   
    % Set up the rw object for this program
    rwObject = IBIOColorDetectReadWriteBasic;
    theProgram = mfilename;

    % Generate simulation params
    [sessionParams, spatialParams, temporalParams, ...
        LMSsamplingParams, chromaticDirectionParams, backgroundParams, ...
        oiParams, mosaicParams, responseSubSamplingParams] = paramsGenerate(instancesNum, temporalResponseToIncludeInUnitsOfIntegrationTime);
    
 
    % Generate the background scene
    [backgroundScene, ~] = generateGaborDisplayScene(spatialParams, backgroundParams, 0.0);
    
    % Generate custom optics
    theOI = colorDetectOpticalImageConstruct(oiParams);
    
    % Compute the background OI
    oiBackground = theOI;
    oiBackground = oiCompute(oiBackground, backgroundScene);
  
    % Generate the coneMosaic or load a previously saved one.
    if (recomputeConeMosaic)
        % Generate the cone mosaic
        theConeMosaic = coneMosaicGenerate(mosaicParams); 
        % Save cone mosaic
        coneParamsList = {mosaicParams};
        rwObject.write('coneMosaic', theConeMosaic, coneParamsList, theProgram, 'type', 'mat');
    else
        % Load the cone mosaic
        coneParamsList = {mosaicParams};
        theConeMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    end
    
    % Retrieve cone indices of interest
    activeConesNum = numel(find(theConeMosaic.pattern > 1));
    [lConeIndices, mConeIndices, sConeIndices, centerMostLconeIndex, centerMostMconeIndex, centerMostSconeIndex] = ...
                    retrieveConeIndices(theConeMosaic);
         
    % Load the 
    conditionIndex = 0;
    for chromaticDirectionIndex = 1:numel(chromaticDirectionParams) 
        % Load params for this chromaticDirection
        chromaticParams =  struct(...
             'backgroundxyY', backgroundParams.backgroundxyY, ...
             'coneContrasts', chromaticDirectionParams{chromaticDirectionIndex}.coneContrastUnitVector);
         
        stimulusStrengthAxis = LMSsamplingParams{chromaticDirectionIndex}.stimulusStrengthAxis;
        instancesNum = LMSsamplingParams{chromaticDirectionIndex}.instancesNum;
          
        fprintf('Condition %d has %d examined strengths\n', conditionIndex, numel(stimulusStrengthAxis));
        
        for stimStrengthIndex = 1:numel(stimulusStrengthAxis)
            % Update the current stimulus strength
            chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrength = stimulusStrengthAxis(stimStrengthIndex);
            
            % Compute modulated stimulus scene
            [modulatedScene, actualStimulusStrength] = generateGaborDisplayScene(spatialParams, chromaticParams, chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrength);
            fprintf('Stimulus %d/%d unit vector: [%2.1f %2.2f %2.2f], Strength (RMS cone contrast): specified = %2.3f, measured: %2.3f\n', ...
                stimStrengthIndex, numel(stimulusStrengthAxis), ...
                chromaticDirectionParams{chromaticDirectionIndex}.coneContrastUnitVector(1), ...
                chromaticDirectionParams{chromaticDirectionIndex}.coneContrastUnitVector(2), ...
                chromaticDirectionParams{chromaticDirectionIndex}.coneContrastUnitVector(3), ...
                stimulusStrengthAxis(stimStrengthIndex), actualStimulusStrength);
            
            % Compute modulated OI
            oiModulated = theOI;
            oiModulated = oiCompute(oiModulated, modulatedScene);
            
            % Generate stimulus label
            stimLabel = sprintf('AzimuthAngle=%2.1fDeg, ElevationAngle=%2.1fDeg, StimStrength=%2.3f', ...
                LMSsamplingParams{chromaticDirectionIndex}.azimuthAngle, ...
                LMSsamplingParams{chromaticDirectionIndex}.elevationAngle, ...
                chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrength);
            
            % Generate response storage struct
            responseData = struct(...
                'stimLabel', stimLabel, ...
                'stimConeContrastUnitVector', chromaticDirectionParams{chromaticDirectionIndex}.coneContrastUnitVector, ...
                'stimStrength',  chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrength, ...
                'absorptionsCountSequence', [], ...
                'absorptionsTimeAxis', [], ... 
                'photoCurrentSignals', [], ... 
                'photoCurrentTimeAxis', []);
            
            % Generate ancillary data struct
            ancillaryData = struct(...
                'sessionParams', sessionParams, ...
                'oiParams', oiParams, ...
                'mosaicParams', mosaicParams, ...
                'spatialParams', spatialParams, ...
                'temporalParams', temporalParams, ...
                'backgroundParams',  backgroundParams, ...
                'LMSsamplingParams',  LMSsamplingParams{chromaticDirectionIndex}, ...
                'chromaticDirectionParams', chromaticDirectionParams{chromaticDirectionIndex}, ...
                'responseSubSamplingParams', responseSubSamplingParams ...
            );
        
            % Params list to define filepath for saved data
            paramsList = {mosaicParams, sessionParams, oiParams, spatialParams, temporalParams, backgroundParams, ...
                          LMSsamplingParams{chromaticDirectionIndex}, chromaticDirectionParams{chromaticDirectionIndex}, responseSubSamplingParams};
                  
            
            dataStruct = struct();
            dataStruct.chromaticDirectionIndex = chromaticDirectionIndex;
            dataStruct.stimStrengthIndex = stimStrengthIndex;
            dataStruct.randomSeed = randi(1000000)+1;
            dataStruct.theConeMosaic = theConeMosaic.copy();   % make a copy because the @coneMosaic is a handle
            dataStruct.theEMpaths = [];
            dataStruct.theOIsequence =  [];
            dataStruct.oiModulated = oiModulated;
            dataStruct.paramsList = paramsList;
            dataStruct.stimLabel = stimLabel;
            dataStruct.responseData = responseData;
            dataStruct.ancillaryData = ancillaryData;
            dataStruct.rwObject = rwObject;

            conditionIndex = conditionIndex + 1;
            parforStructData{conditionIndex} = dataStruct;
        end % stimStrengthIndex
    end % chromaticDirectionIndex
    
    parforConditionsNum = numel(parforStructData); 

    fprintf('Staring parfor loop\n');
    
    % Loop over all conditions
    parfor parforCondIndex = 1:parforConditionsNum
        
        t = getCurrentTask();
        fprintf('<strong>Processor [%d] working on condition %d / %d\n</strong>', t.ID, parformCondIndex, parforConditionsNum);
        
        % Get the data for the current condition
        theData = parforStructData{parforCondIndex};
        
        % Print all params
        if (paramsVerbosity > 0)
            disp(UnitTest.displayNicelyFormattedStruct(theData.ancillaryData, 'ancillaryData', '', 60));
        end
            
        % Compute the oiSequence
        theData.theOIsequence = oiSequence(oiBackground, theData.oiModulated, theData.ancillaryData.temporalParams.stimTimeAxis, ...
                                   theData.ancillaryData.temporalParams.stimulusModulationFunction, 'composition', 'blend');
        %theData.theOIsequence.visualize('format', 'montage');
        
        % Initialize random number generator
        rng(theData.randomSeed);
        
        % Compute the emPaths for this condition
        eyeMovementsNo = computeEyeMovementsNum(theData.theConeMosaic.integrationTime, theData.theOIsequence);

        for instanceIndex = 1:instancesNum
            theData.theEMpaths(instanceIndex, :,:) = theData.theConeMosaic.emGenSequence(eyeMovementsNo);
        end
        if (theData.ancillaryData.mosaicParams.eyesDoNotMove)
            theData.theEMpaths(instanceIndex, :,:) = 0*theData.theEMpaths(instanceIndex, :,:);
        end
           
        % Compute the responses
        fprintf('%s\nComputing %d response instances for a %d x %d cone mosaic (active cones:%d) with %d eye movements scanning a %d-frame stimulus sequence.\n', ...
                theData.stimLabel, instancesNum, theData.theConeMosaic.mosaicSize(1), theData.theConeMosaic.mosaicSize(2), activeConesNum, eyeMovementsNo, theData.theOIsequence.length);
           
        
        tic
        responseData = theData.responseData;
        
        % Compute all responses instances for this condition
        [responseData.absorptionsCountSequence, ...
         responseData.absorptionsTimeAxis, ...
         responseData.photoCurrentSignals, ...
         responseData.photoCurrentTimeAxis] = ...
                    theData.theConeMosaic.computeForOISequence(theData.theOIsequence, ...
                    'emPaths', theData.theEMpaths, ...
                    'currentFlag', true, ...
                    'newNoise', true ...
                    );
        fprintf('Response computation took %2.2f minutes\n', toc/60);
        
        % Determine portion of the absorptions signal to be kept
        absorptionsTimeIndicesToKeep = determineTimeIndicesToKeep(responseData.absorptionsTimeAxis, ...
            theData.theConeMosaic.integrationTime, theData.ancillaryData.temporalParams.rampPeak, ...
            theData.ancillaryData.responseSubSamplingParams);
        
        % Determine portion of the photocurrents signal to be kept
        photocurrentsTimeIndicesToKeep = determineTimeIndicesToKeep(responseData.photoCurrentTimeAxis, ...
            theData.theConeMosaic.integrationTime, theData.ancillaryData.temporalParams.rampPeak, ...
            theData.ancillaryData.responseSubSamplingParams);
        
        % Remove unwanted portions of the time axes
        responseData.absorptionsTimeAxis = responseData.absorptionsTimeAxis(absorptionsTimeIndicesToKeep);
        responseData.photoCurrentTimeAxis = responseData.photoCurrentTimeAxis(photocurrentsTimeIndicesToKeep);
        
        % Remove unwanted portions of the responses
        if (isa(theData.theConeMosaic, 'coneMosaicHex'))
            responseData.absorptionsCountSequence = responseData.absorptionsCountSequence(:,:,absorptionsTimeIndicesToKeep);
            responseData.photoCurrentSignals = responseData.photoCurrentSignals(:,:,photocurrentsTimeIndicesToKeep);
        else
            responseData.absorptionsCountSequence = responseData.absorptionsCountSequence(:,:,:,absorptionsTimeIndicesToKeep);
            responseData.photoCurrentSignals = responseData.photoCurrentSignals(:,:,:,photocurrentsTimeIndicesToKeep);
        end
        
        % Save responses and parameters for this stimStrengthIndex and chromaticDirectionIndex
        theData.rwObject.write('coneResponses', responseData, theData.paramsList, theProgram, 'type', 'mat');
            
        % Save other data we need for use by the classifier preprocessing subroutine
        theData.rwObject.write('ancillaryData', theData.ancillaryData, theData.paramsList, theProgram, 'type', 'mat');
            
        % Export stills and video of 2D mosaic activation (absorptions + photocurrents) for some instances
        if (exportMosaic2DActivationStillsAndVideo)
            if (strcmp(theData.ancillaryData.mosaicParams.mosaicType, 'HEX')) && ~any(isnan(theData.ancillaryData.mosaicParams.fieldOfViewDegs))  
                instancesToVisualize = 10;
                exportHexMosaicActivationStillsAndVideos(responseData, theData.theConeMosaic, instancesToVisualize, theData.rwObject, theData.paramsList, theProgram);
            end
        end
            
        % Export stills of time traces (absorptions + photocurrents)
        if (exportLMSresponseTraceStills)
            coneStride = 50;
            exportResponseTraceStills(theData.theConeMosaic, reponseData, ...
                theData.ancillaryData.chromaticDirectionParams{chromaticDirectionIndex}, theData.ancillaryData.temporalParams, coneStride, ...
                lConeIndices, mConeIndices, sConeIndices, centerMostLconeIndex, centerMostMconeIndex, centerMostSconeIndex, ...
                theData.rwObject, theData.paramsList, theProgram);
        end 
    end % all conditions parfor
end


% ------- Helper functions ------

function [sessionParams, spatialParams, temporalParams, ...
        LMSsamplingParams, chromaticDirectionParams, backgroundParams, ...
        oiParams, mosaicParams, responseSubSamplingParams] = paramsGenerate(temporalResponseToIncludeInUnitsOfIntegrationTime, instancesNum)
    
    sessionParams = struct(...
                 'type', 'Session', ...
          'sessionType', 'LMSplane' ...
    );

    % In P&W 1996, in the constant cycle condition, this was 10 deg (Section 2.2, p 517)
    PW96_fovDegs = 10;
    
    spatialParams = struct(...
                  'type', 'Spatial_v2', ...
            'spatialType', 'Gabor', ...
        'fieldOfViewDegs', PW96_fovDegs, ... 
        'cyclesPerDegree', 2.0,...
      'spatialPhaseDegs', 0, ...
       'orientationDegs', 90, ...
             'windowType', 'Gaussian', ...
        'gaussianFWHMDegs', 1.9, ...
        'viewingDistance', 0.75, ...
                    'row', 256, ...
                    'col', 256 ...
        );

    
    CRTrefreshRate = 87;
    temporalParams = struct(...
                         'type', 'Temporal_v2', ...
                    'frameRate', CRTrefreshRate, ...
 'stimulusSamplingIntervalSecs', 1.0/CRTrefreshRate, ...
             'rampDurationSecs', 900/1000, ... 
                  'rampTauSecs', 165/1000, ...
                     'rampPeak', 0/1000, ... 
      'stimTimeAxisOriginShift', -100/1000 ...   % allow this many milliseconds for the photocurrent response to stabilize
        );
    
    
    % Define the stimulus temporal modulation function
    temporalParams.stimTimeAxis = 0:temporalParams.stimulusSamplingIntervalSecs:temporalParams.rampDurationSecs;
    temporalParams.stimTimeAxis = temporalParams.stimTimeAxisOriginShift + temporalParams.stimTimeAxis - mean(temporalParams.stimTimeAxis);
    % Compute modulation function
    temporalParams.stimulusModulationFunction = exp(-0.5*((temporalParams.stimTimeAxis-temporalParams.rampPeak)/temporalParams.rampTauSecs).^2);
    
    % Optics params
    oiParams = oiParamsGenerate();
    oiParams.type = 'Optics_v2';
    oiParams.fieldOfViewDegs = spatialParams.fieldOfViewDegs*1.2;
    
    
    mosaicParams = struct(...
                'type', 'Mosaic_v2', ...
          'mosaicType', 'HEX', ...
     'fieldOfViewDegs', PW96_fovDegs*[0.3 0.3],...         % nan for 1L, 1M, and 1S-cone only
    'eccentricityDegs', 0, ...  
 'spatialLMSDensities', [0.62 0.31 0.07], ...
 'integrationTimeSecs', 50/1000, ...            % 50 msec integration time
         'photonNoise', true, ...              % add Poisson noise
      'osTimeStepSecs', 10/1000, ...             % 5 milliseconds
             'osNoise', true, ...              % outer-segment noise
       'eyesDoNotMove', false ....              % normal eye movements
       );
   
   % Response subSampling (how many seconds to include)
   responseSubSamplingParams = struct(...
                     'type', 'ResponseSubsampling', ...
         'secondsToInclude', temporalResponseToIncludeInUnitsOfIntegrationTime*mosaicParams.integrationTimeSecs, ...          % temporal subsampling: only keep responses within this time period around 
   'secondsToIncludeOffset', mosaicParams.integrationTimeSecs*0.5 ...   % 'temporalParams.rampPeak', and this offset
        );
   
    % In the constant cycle condition, the background was xyY= 0.38, 0.39, 536.2 cd/m2
    % Also they say they placed a uniform field to the screen to increase the contrast resolutions (page 517, Experimental Aparratus section)
    xyY =  [0.38 0.39 536.2];
    backgroundParams = struct(...
                     'type', 'Background_v2', ...
                   'device', 'Monitor', ...
            'backgroundxyY', xyY, ...
            'coneContrasts', [0.0 0.0 0.0] ...
        );
    
    % Define we sample sensitivity along different directions in the LMS space  
    LMSsamplingParams = LMSsamplingParamsGenerate(instancesNum);
 
    % Generate chromaticDirectionParams
    chromaticDirectionParams = chromaticDirectionParamsGenerate(LMSsamplingParams);
end
    
    
function LMSsamplingParams = LMSsamplingParamsGenerate(instancesNum)
    % Add sampling params for each of the LMS directions we want to explore
    LMSsamplingParams = {};
    
    % the null stimulus (zero stimulus strength)
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 0, ...                      % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                    % z-axis (S-modulation)
    'stimulusStrengthAxis', 0, ...      
            'instancesNum', instancesNum ...
        );
    
    % L direction, specifically: cL = 1, cM = 0, cS = 0.0
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 0.0, ...                        % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                        % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.01, 0.1, 10), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
            'instancesNum', instancesNum ...
        );
    
    % L+M direction, specifically: cL = 0.7071, cM = 0.7071, cS = 0.0
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 45.0, ...                        % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                        % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.01, 0.1, 10), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
            'instancesNum', instancesNum ...
        );
    
    % M direction, specifically: cL = 0, cM = 1.0, cS = 0.0
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 90.0, ...                        % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                         % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.01, 0.1, 10), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
            'instancesNum', instancesNum ...
        );
    
    
    % -L+M direction, specifically: cL = -0.7071, cM = 0.7071, cS = 0.0
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 135.0, ...                      % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                        % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.01, 0.1, 10), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
            'instancesNum', instancesNum ...
        );  
    
    % S direction, specifically: cL = -0.0, cM = 0.0, cS = 1.0
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 0.0, ...                         % (x/y plane) (L/M modulation)
          'elevationAngle', 90.0, ...                        % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.01, 0.7, 10), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
            'instancesNum', instancesNum ...
        );  
end
    

function indicesToKeep = determineTimeIndicesToKeep(timeAxis, integrationTime, rampPeak, responseSubSamplingParams)
    t = timeAxis + integrationTime/2 - rampPeak - responseSubSamplingParams.secondsToIncludeOffset;
    indicesToKeep  = find(abs(t) <= theData.ancillaryData.responseSubSamplingParams.secondsToInclude/2);
    if (isempty(indicesToKeep ))
       [~,indicesToKeep ] = min(abs(t-responseSubSamplingParams.secondsToInclude/2));
    end
end
        
function exportResponseTraceStills(theConeMosaic, stimData, chromaticDirectionParams, temporalParams, coneStride, lConeIndices, mConeIndices, sConeIndices, centerMostLconeIndex, centerMostMconeIndex, centerMostSconeIndex, rwObject, paramsList, theProgram)
            
    % Plotting time limits
    timeLimits = [stimData.absorptionsTimeAxis(1) stimData.absorptionsTimeAxis(end)];
    if (timeLimits(2) == timeLimits(1))
        timeLimits = timeLimits(1) + theConeMosaic.integrationTime*[-2 2];
    end

    % Plot all absorption response instances for the center-most L-, M-, and S-cone
    instancesNum = size(stimData.absorptionsCountSequence, 1);
    instancesToPlot = 1:instancesNum;
    hFig = plotResponseTimeSeries(...
      'absorptions', ...
      stimData.absorptionsTimeAxis, ...
      stimData.absorptionsCountSequence(instancesToPlot,:,:,:), ...
      theConeMosaic.integrationTime, ...
      temporalParams.stimTimeAxis, timeLimits, temporalParams.stimulusModulationFunction, ...
      centerMostLconeIndex, centerMostMconeIndex, centerMostSconeIndex, 1e12, ...
      chromaticDirectionParams.coneContrastUnitVector, ...
      chromaticDirectionParams.stimulusStrength, []);

    % Plot all photocurrent response instances for the center-most L,M, and S-cone
    plotResponseTimeSeries(...
      'photocurrents', ...
      stimData.photoCurrentTimeAxis, ...
      stimData.photoCurrentSignals(instancesToPlot,:,:,:), ...
      theConeMosaic.integrationTime, ...
      temporalParams.stimTimeAxis, timeLimits, temporalParams.stimulusModulationFunction, ...
      centerMostLconeIndex, centerMostMconeIndex, centerMostSconeIndex, 1e12, ...
      chromaticDirectionParams.coneContrastUnitVector, ...
      chromaticDirectionParams.stimulusStrength, hFig);

    % Export a PNG image    
    rwObject.write('AbsorptionsPhotocurrentsSingleConesAllInstances', stimData, paramsList, theProgram, ...
        'type', 'NicePlotExport', 'FigureHandle', hFig, 'FigureType', 'png');

    % Plot some cone responses (number of cones is controlled by the coneStride parameter) for the first instance
    if (1==1)
        % Plot some L,M, and S cone absorption responses for the first few intances 
        instancesToPlot = 1:1;

        hFig = plotResponseTimeSeries(...
          'absorptions', ...
          stimData.absorptionsTimeAxis, ...
          stimData.absorptionsCountSequence(instancesToPlot,:,:,:), ...
          theConeMosaic.integrationTime, ...
          temporalParams.stimTimeAxis, timeLimits, temporalParams.stimulusModulationFunction, ...
          lConeIndices, mConeIndices, sConeIndices, coneStride, ...
          chromaticDirectionParams.coneContrastUnitVector, ...
          chromaticDirectionParams.stimulusStrength, []);

       % Plot some L,M, and S cone photocurrent responses for the first few intances
        plotResponseTimeSeries(...
          'photocurrents', ...
          stimData.photoCurrentTimeAxis, ...
          stimData.photoCurrentSignals(instancesToPlot,:,:,:), ...
          theConeMosaic.integrationTime, ...
          temporalParams.stimTimeAxis, timeLimits, temporalParams.stimulusModulationFunction, ...
          lConeIndices, mConeIndices, sConeIndices,  coneStride, ...
          chromaticDirectionParams.coneContrastUnitVector, ...
          chromaticDirectionParams.stimulusStrength, hFig);

     % Export a PNG image 
     rwObject.write('AbsorptionsPhotocurrentsSelectConesSelectInstances', stimData, paramsList, theProgram, ...
         'type', 'NicePlotExport', 'FigureHandle', hFig, 'FigureType', 'png');
    end
end
            
function [iL, iM, iS, lconeToPlot, mconeToPlot, sconeToPlot]  = retrieveConeIndices(theConeMosaic)
    % Find L, M, and S-cone indices for response plotting
    iL = find(theConeMosaic.pattern == 2);
    iM = find(theConeMosaic.pattern == 3);
    iS = find(theConeMosaic.pattern == 4);

    % Find center-most L cone (row,col) coord
    [~,idxL] = min(sum(theConeMosaic.coneLocs(iL,:).^2,2));
    % Find center-most M cone (row,col) coord
    [~,idxM] = min(sum(theConeMosaic.coneLocs(iM,:).^2,2));
    % Find center-most S cone (row,col) coord
    [~,idxS] = min(sum(theConeMosaic.coneLocs(iS,:).^2,2));

    % Find indices for responses for the center-most L,M and S-cone
    % This works also with a hex mosaic, which returns 
    % responses for only the non-null cones
    nonNullConeIndices = find(theConeMosaic.pattern > 1);
    lconeToPlot = find(nonNullConeIndices == iL(idxL));
    mconeToPlot = find(nonNullConeIndices == iM(idxM));
    sconeToPlot = find(nonNullConeIndices == iS(idxS));

    % Find indices for L,M,S cone responses
    % This works also with a hex mosaic, which returns 
    % responses for only the non-null cones
    iL = find(theConeMosaic.pattern(nonNullConeIndices) == 2);
    iM = find(theConeMosaic.pattern(nonNullConeIndices) == 3);
    iS = find(theConeMosaic.pattern(nonNullConeIndices) == 4);
end
                
function exportHexMosaicActivationStillsAndVideos(stimData, theConeMosaic, instancesToVisualize,  rwObject, paramsList, theProgram)
            
    videoAbsorptionsFilename = 'absorptionsVideo';
    videoAbsorptionsOBJ = VideoWriter(videoAbsorptionsFilename, 'MPEG-4'); % H264 format
    videoAbsorptionsOBJ.FrameRate = 30; 
    videoAbsorptionsOBJ.Quality = 100;
    videoAbsorptionsOBJ.open();

    videoPhotocurrentsFilename = 'photoCurrentsVideo';
    videoPhotocurrentsOBJ = VideoWriter(videoPhotocurrentsFilename, 'MPEG-4'); % H264 format
    videoPhotocurrentsOBJ.FrameRate = 30; 
    videoPhotocurrentsOBJ.Quality = 100;
    videoPhotocurrentsOBJ.open();

    % Determine signal ranges
    absorptionsRange = [min(stimData.absorptionsCountSequence(:)) max(stimData.absorptionsCountSequence(:))];
    photocurrentsRange = [min(stimData.photoCurrentSignals(:)) max(stimData.photoCurrentSignals(:))];   

    aspectRatio = 1.0;
    % Initial zoom-in
    zoomInFactorAbsorptions = aspectRatio;
    zoomInFactorPhotocurrents = aspectRatio;

    activationLUT = bone(1024);

    for visualizedInstanceIndex = 1:min([size(stimData.absorptionsCountSequence,1) instancesToVisualize])

        visualizedAbsorptions = squeeze(stimData.absorptionsCountSequence(visualizedInstanceIndex,:,:,:));
        sz = ndims(visualizedAbsorptions);
        absorptionSequenceLength = size(visualizedAbsorptions, sz(end));
        
        for tIndex = 1:absorptionSequenceLength
            % Compute zoom-in factor
            if (visualizedInstanceIndex>3)
                zoomInFactorAbsorptions = zoomInFactorAbsorptions * (1.0 - 0.05/absorptionSequenceLength);
            end

            if (zoomInFactorAbsorptions < 0.4)
                zoomInFactorAbsorptions = 0.4;
            end

            if (tIndex > 1) 
                close(hFig);
            end

            hFig = theConeMosaic.visualizeActivationMaps(...
                squeeze(visualizedAbsorptions(:,tIndex)), ...                                           % the signal matrix
                   'mapType', 'modulated disks', ...                                    % how to display cones: choose between 'density plot', 'modulated disks' and 'modulated hexagons'
                'signalName', 'isomerizations (R*/cone/integration time)', ...          % colormap title (signal name and units)
               'signalRange', absorptionsRange, ...                                    % signal range
                    'xRange', theConeMosaic.width/2 * [-1 1], ...
                    'yRange', theConeMosaic.height/2 * [-1 1]*0.78, ...
              'zoomInFactor', zoomInFactorAbsorptions, ...                              % dynamically changing zoom
                  'colorMap', activationLUT, ...                                        % colormap to use for displaying activation level
        'separateLMSmosaics', false, ...                                                % when true, L-,M-, and S-cone submosaic activations are plotted separately
            'activationTime', stimData.absorptionsTimeAxis(tIndex), ...                 % current response time
   'visualizedInstanceIndex', visualizedInstanceIndex, ...                              % currently visualized instance
                'figureSize', [1280 960] ...                                            % figure size in pixels
            );

            % Export a PNG image    
            filename = sprintf('MosaicAbsorptions_%2.3f', stimData.absorptionsTimeAxis(tIndex));
            rwObject.write(filename, stimData, paramsList, theProgram, ...
                    'type', 'NicePlotExport', 'FigureHandle', hFig, 'FigureType', 'png');

            videoAbsorptionsOBJ.writeVideo(getframe(hFig));
        end % tIndex
        close(hFig);
        
        videoAbsorptionsOBJ.close();
        % Export video to right directory
        rwObject.write(videoAbsorptionsFilename, fullfile(pwd,sprintf('%s.mp4',videoAbsorptionsFilename)), ...
        paramsList,theProgram,'Type','movieFile', 'MovieType', 'mp4');
    

        % Photocurrents video
        visualizedPhotocurrents = squeeze(stimData.photoCurrentSignals(visualizedInstanceIndex,:,:,:));
        sz = ndims(visualizedPhotocurrents);
        photocurrentsSignalLength = size(visualizedPhotocurrents, sz(end));
        
        for tIndex = 1:photocurrentsSignalLength
            % compute zoom-in factor
            if (visualizedInstanceIndex > 3)
                zoomInFactorPhotocurrents = zoomInFactorPhotocurrents * (1.0 - 0.05/photocurrentsSignalLength);
            end

            if (zoomInFactorPhotocurrents < 0.4)
                zoomInFactorPhotocurrents = 0.4;
            end
            if (tIndex > 1) 
                close(hFig2);
            end
            hFig2 = theConeMosaic.visualizeActivationMaps(...
                squeeze(visualizedPhotocurrents(:,tIndex)), ...                                           % the signal matrix
                   'mapType', 'modulated disks', ...                        % how to display cones: choose between 'density plot', 'modulated disks' and 'modulated hexagons'
                'signalName', 'photocurrent (pAmps)', ...                   % colormap title (signal name and units)
               'signalRange', photocurrentsRange, ...                       % signal range
                    'xRange', theConeMosaic.width/2 * [-1 1], ...
                    'yRange', theConeMosaic.height/2 * [-1 1]*0.78, ...
              'zoomInFactor', zoomInFactorPhotocurrents, ...                % dynamically changing zoom
                  'colorMap', activationLUT, ...                            % colormap to use for displaying activation level
        'separateLMSmosaics', false, ...                                    % when true, L-,M-, and S-cone submosaic activations are plotted separately
            'activationTime', stimData.photoCurrentTimeAxis(tIndex), ...    % current response time
   'visualizedInstanceIndex', visualizedInstanceIndex, ...                  % currently visualized instance
                'figureSize', [1280 960] ...                                % figure size in pixels
            );

            % Export a PNG image    
            filename = sprintf('MosaicPhotocurrents_%2.3f', stimData.photoCurrentTimeAxis(tIndex));

            rwObject.write(filename, stimData, paramsList, theProgram, ...
                    'type', 'NicePlotExport', 'FigureHandle', hFig2, 'FigureType', 'png');

            videoPhotocurrentsOBJ.writeVideo(getframe(hFig2));
        end % tIndex
        
    end % visualizedInstanceIndex
    close(hFig2);
    
    videoPhotocurrentsOBJ.close();
    % Export video to right directory 
    rwObject.write(videoPhotocurrentsFilename, fullfile(pwd,sprintf('%s.mp4',videoPhotocurrentsFilename)), ...
        paramsList,theProgram,'Type','movieFile', 'MovieType', 'mp4');
end
                
function hFig = plotResponseTimeSeries(signalName, timeAxis, responseTimeSeries, integrationTime, stimTimeAxis, timeLimits, ...
                stimulusModulationFunction, iL, iM, iS,  coneStride, coneContrastUnitVector, stimulusStrength, hFig0)
            
    instancesPlotted = size(responseTimeSeries,1);
    timePointsPlotted = numel(timeAxis);
    
    minResponseTimeSeries = min(responseTimeSeries(:));
    maxResponseTimeSeries = max(responseTimeSeries(:));
    
    if (strcmp(signalName, 'photocurrents'))
        minResponseTimeSeries = -70;
        maxResponseTimeSeries = 0;
    end
    
    % Extract cone responses by type
    iL = iL(1:coneStride:end);
    iM = iM(1:coneStride:end);
    iS = iS(1:coneStride:end);
    
    if (timePointsPlotted == 1)
        responseTimeSeriesByConeType{1} = reshape(responseTimeSeries(:,iL), [instancesPlotted*numel(iL) 1]);
        responseTimeSeriesByConeType{2} = reshape(responseTimeSeries(:,iM), [instancesPlotted*numel(iM) 1]);
        responseTimeSeriesByConeType{3} = reshape(responseTimeSeries(:,iS), [instancesPlotted*numel(iS) 1]);
    else
        if (ndims(responseTimeSeries) == 4)
            responseTimeSeries = reshape(responseTimeSeries, [instancesPlotted size(responseTimeSeries,2)*size(responseTimeSeries,3) timePointsPlotted]);
        end
        responseTimeSeries = permute(responseTimeSeries, [1 3 2]);
        responseTimeSeriesByConeType{1} = reshape(permute(responseTimeSeries(:,:,iL), [1 3 2]), [instancesPlotted*numel(iL) timePointsPlotted]);
        responseTimeSeriesByConeType{2} = reshape(permute(responseTimeSeries(:,:,iM), [1 3 2]), [instancesPlotted*numel(iM) timePointsPlotted]);
        responseTimeSeriesByConeType{3} = reshape(permute(responseTimeSeries(:,:,iS), [1 3 2]), [instancesPlotted*numel(iS) timePointsPlotted]);
    end
    colors = [1 0 0; 0 1.0 0; 0 0.8 1];
    plotBackgroundColor = [0.1 0.1 0.1];
    barOpacity = 0.2;
    
    if (isempty(hFig0))
        hFig = figure(); clf;
        set(hFig, 'Color', [0 0 0], 'Position', [10 10 1500 900]);
        subplot('Position', [0.06 0.06 0.42 0.87]);
    else
       figure(hFig0); 
       hFig = hFig0;
       subplot('Position', [0.54 0.06 0.42 0.87]);
    end
    
    yyaxis left
    hold on
    % Identify stimulus presentation times
    for k = 1:numel(stimTimeAxis)
        plot(stimTimeAxis(k)*[1 1]*1000, [minResponseTimeSeries maxResponseTimeSeries], 'k-', 'Color', [0.5 0.5 0.5]);
    end
        
    dt = integrationTime;
    % Plot responseTimeSeries
    switch signalName
        case 'absorptions'
            for coneType = 1:3
                for tIndex = 1:numel(timeAxis)
                    quantaAtThisTimeBin = squeeze(responseTimeSeriesByConeType{coneType}(:,tIndex));
                    plot([timeAxis(tIndex) timeAxis(tIndex)+dt]*1000, [quantaAtThisTimeBin(:) quantaAtThisTimeBin(:)], '-', ...
                        'LineWidth', 2.0, 'Color', [colors(coneType,:) barOpacity]);
                end
            end % coneType
            
        case 'photocurrents'
            for coneType = 1:3
                signalsPlotted = size(responseTimeSeriesByConeType{coneType},1);
                for signalIndex = 1:signalsPlotted
                    plot(timeAxis*1000, squeeze(responseTimeSeriesByConeType{coneType}(signalIndex,:)), '-', 'LineWidth', 2.0, 'Color', [colors(coneType,:) barOpacity]);
                end
            end % coneType
        otherwise
            error('Unknown signal name: ''%s''.', signalName);
    end
    
    box on;
    xlabel('time (ms)', 'FontSize', 14);
    switch signalName
        case 'absorptions'
            ylabel(sprintf('# of absorptions per %2.1f ms',dt*1000), 'FontSize', 14);
        case 'photocurrents'
            ylabel(sprintf('photocurrent (pA)'), 'FontSize', 14);
        otherwise
            error('Unknown signal name: ''%s''.', signalName);
    end
    
    set(gca, 'XLim', timeLimits*1000, ...
             'YLim', [minResponseTimeSeries maxResponseTimeSeries], ...
             'Color', plotBackgroundColor, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], ...
             'FontSize', 14);
       
    % Plot the stimulus modulation on the right time axis
    yyaxis right
    dt = stimTimeAxis(2)-stimTimeAxis(1);
    for tIndex = 1:numel(stimTimeAxis)
        plot([stimTimeAxis(tIndex) stimTimeAxis(tIndex)+dt]*1000, stimulusModulationFunction(tIndex)*[1 1], '-', 'LineWidth', 1.5, 'Color', 'y');
        if (tIndex < numel(stimTimeAxis))
            plot(stimTimeAxis(tIndex+1)*[1 1]*1000, [stimulusModulationFunction(tIndex) stimulusModulationFunction(tIndex+1)], '-', 'LineWidth', 1.5, 'Color', 'y');
        end
    end
    set(gca, 'YColor', 'y', 'FontSize', 14);
    if (~isempty(hFig0))
        ylabel('stimulus modulation', 'FontSize', 14);
    end

    title(sprintf('LMS unit vector: <%0.2f, %0.2f, %0.2f>, stim. strength: %0.5f\n %d response instances from %d L-, %d M-, and %d S-cone(s)', ...
        coneContrastUnitVector(1), coneContrastUnitVector(2), coneContrastUnitVector(3), stimulusStrength, ...
        instancesPlotted, numel(iL), numel(iM), numel(iS)), 'Color', [1 1 1], 'FontSize',16, 'FontWeight', 'normal');
    drawnow;
end


function chromaticDirectionParams = chromaticDirectionParamsGenerate(LMSsamplingParams)
    
    % Initialize the chromaticDirectionParams cell array
    chromaticDirectionParams = cell(1,numel(LMSsamplingParams));
    
    % Populate the chromaticDirectionParams
    for chromaticDirectionIndex = 1:numel(LMSsamplingParams)  
        % Get angle on LM modulation plane
        azimuthAngle = LMSsamplingParams{chromaticDirectionIndex}.azimuthAngle;
        % Get angle along S-modulation axis
        elevationAngle = LMSsamplingParams{chromaticDirectionIndex}.elevationAngle;
            
         % Compute cone contrasts from azimuth and elevation
        [cL, cM, cS] = sph2cart(azimuthAngle/180*pi, elevationAngle/180*pi, 1);

        % Normalize to unity RMS cone contrast (stimulus strength)
        coneContrasts = [cL, cM, cS];
        coneContrastUnitVector = coneContrasts / norm(coneContrasts);

        % Form chromaticDirection params
        chromaticDirectionParams{chromaticDirectionIndex} = struct(...
                          'type', 'ColorModulation_v2', ...
                        'device', 'Monitor', ...
        'coneContrastUnitVector', coneContrastUnitVector, ...
              'stimulusStrength', [] ...           % the current stimulus strength,  updated at runtime
                );
    end  % chromaticDirectionIndex
end

    
function theConeMosaic = coneMosaicGenerate(mosaicParams)

    if (strcmp(mosaicParams.mosaicType, 'HEX')) && ~any(isnan(mosaicParams.fieldOfViewDegs))
        % HEX mosaic
        resamplingFactor = 6;
        centerInMM = [0.0 0.0];                     % mosaic eccentricity
        spatiallyVaryingConeDensity = true;        % constant spatial density (at the mosaic's eccentricity)

        theConeMosaic = coneMosaicHex(resamplingFactor, spatiallyVaryingConeDensity, ...
                       'center', centerInMM*1e-3, ...
               'spatialDensity', [0 mosaicParams.spatialLMSDensities] ...
            );
        theConeMosaic.setSizeToFOVForHexMosaic(mosaicParams.fieldOfViewDegs);
        theConeMosaic.visualizeGrid();
        
    else
        % RECT mosaic
        theConeMosaic = coneMosaic;

        % Adjust size
        if isnan(mosaicParams.fieldOfViewDegs)
            % Generate a human cone mosaic with 1L, 1M and 1S cone
            theConeMosaic.rows = 1;
            theConeMosaic.cols = 3;
            theConeMosaic.pattern = [2 3 4];
        else
            theConeMosaic.setSizeToFOV(mosaicParams.fieldOfViewDegs);
            % Set the LMS spatial densities
            theConeMosaic.spatialDensity = [0 mosaicParams.spatialLMSDensities]';
        end
    end

    % Set the noise
    theConeMosaic.noiseFlag = mosaicParams.photonNoise;

    % Set the integrationTime
    theConeMosaic.integrationTime = mosaicParams.integrationTimeSecs;
    
    % Generate the outer-segment object to be used by the coneMosaic
    theOuterSegment = osLinear();
    theOuterSegment.noiseFlag = mosaicParams.osNoise;
    
    % Set a custom timeStep, for @osLinear we do not need the default 0.1 msec
    theOuterSegment.timeStep = mosaicParams.osTimeStepSecs;

    % Couple the outersegment object to the cone mosaic object
    theConeMosaic.os = theOuterSegment;
end


function eyeMovementsNum = computeEyeMovementsNum(integrationTime, theOIsequence)
    % Generate eye movement sequence for all oi's
    stimulusSamplingInterval = theOIsequence.oiTimeAxis(2)-theOIsequence.oiTimeAxis(1);
    eyeMovementsNumPerOpticalImage = stimulusSamplingInterval/integrationTime;
    eyeMovementsNum = round(eyeMovementsNumPerOpticalImage*theOIsequence.length);
    
    if (eyeMovementsNum < 1)
        error('Less than 1 eye movement!!! \nStimulus sampling interval:%g ms Cone mosaic integration time: %g ms\n', 1000*stimulusSamplingInterval, 1000*theConeMosaic.integrationTime);
    else 
        fprintf('Optical image sequence contains %2.0f eye movements (%2.2f eye movements/oi)\n', eyeMovementsNum, eyeMovementsNumPerOpticalImage);
    end 
end

function [displayScene, measuredStimulusStrength] = generateGaborDisplayScene(spatialParams, colorParams, nominalStimulusStrength)

    % Genereate the rendering display
    display = displayCreate('CRT-MODEL');
    display = displaySet(display,'viewingdistance', 1.0);
    
    % Adjust display's SPDs so as to be able to generate the desired luminance
    peakLuminanceBeforeAdjustment = displayGet(display, 'peak luminance');
    spdMultiplierToGetDesiredLum = colorParams.backgroundxyY(3)/(0.5*peakLuminanceBeforeAdjustment);
    display = displaySet(display,'spd',spdMultiplierToGetDesiredLum*displayGet(display,'spd'));
    
    % Get the SPDs and their spectral sampling
    displayWls = displayGet(display,'wave');
    displaySpectralSampling = WlsToS(displayWls);
    displaySPD = displayGet(display,'spd');
    
    % Load cone fundamentals and XYZtoConeExcitationsMatrix
    [XYZtoConeExcitationsMatrix, coneFundamentals, coneSpectralSampling] = XYZToCones();
    
    % Resample cone fundamentals to the spectral sampling of the display
    coneFundamentals = SplineCmf(coneSpectralSampling, coneFundamentals, displaySpectralSampling);
    
    % Compute coneExcitationsToDisplayPrimaryMatrix
    displayPrimaryToConeExcitations = coneFundamentals * displaySPD * displaySpectralSampling(2);
    coneExcitationsToDisplayPrimaryMatrix = inv(displayPrimaryToConeExcitations);
    
    % Convert the background xyY to cone excitations
    backgroundConeExcitations = XYZtoConeExcitationsMatrix * xyYToXYZ(colorParams.backgroundxyY');
    
    % Generate a Gabor spatial modulation pattern
    gaborModulationPattern = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, nominalStimulusStrength))-1;
    
    % Compute cone excitationImage using the gaborModulationPattern, the
    % backgroundConeExcitations and the coneContrasts
    for ii = 1:3
        coneExcitation = backgroundConeExcitations(ii) * (1 + gaborModulationPattern * colorParams.coneContrasts(ii));
        minExc = min(coneExcitation(:));
        maxExc = max(coneExcitation(:));
        actualConeContrast(ii) = (maxExc-minExc)/(maxExc+minExc);
        coneExcitationImage(:,:,ii) = coneExcitation;
    end
    measuredStimulusStrength = sqrt(sum(actualConeContrast.^2));
    
    % Compute the RGB primary values
    [coneExcitations,m,n] = ImageToCalFormat(coneExcitationImage);
    displayRGBPrimaries = coneExcitationsToDisplayPrimaryMatrix * coneExcitations;
    displayRGBPrimaryImage = CalFormatToImage(displayRGBPrimaries,m,n);
    
    % Check for out-of-gamut
    if (any(displayRGBPrimaryImage(:) > 1))
        fprintf(2,'Image above gamut\n');
    end
    if (any(displayRGBPrimaryImage(:) < 0))
        fprintf(2,'Image below gamut\n');
    end
    
    % Gamma correct the primary values
    displayRGBsettingsImage = round(ieLUTLinear(displayRGBPrimaryImage,displayGet(display,'inverse gamma')));
    
    % Make a scene from the gaborRGBsettings values
    displayScene = sceneFromFile(displayRGBsettingsImage,'rgb',[],display);
    displayScene = sceneSet(displayScene, 'h fov', spatialParams.fieldOfViewDegs);
    displayScene = sceneSet(displayScene, 'distance', spatialParams.viewingDistance);
end

function [M, coneFundamentals, coneSpectralSampling] = XYZToCones
    % Here we'll use the Stockman-Sharpe 2-degree fundamentals and the proposed CIE corresponding XYZ functions
    theCones = load('T_cones_ss2');
    theXYZ = load('T_xyzCIEPhys2');
    
    XYZcolorMatchingFunctions = 683 * theXYZ.T_xyzCIEPhys2;
    coneFundamentals = 683 * theCones.T_cones_ss2;
    coneSpectralSampling = theCones.S_cones_ss2;
    M = ((XYZcolorMatchingFunctions')\(coneFundamentals'))';
end