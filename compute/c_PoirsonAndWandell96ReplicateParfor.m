function c_PoirsonAndWandell96ReplicateParfor
% c_PoirsonAndWandell96Replicate
%
% Compute color detection thresholds to replicate the Poirson & Wandell 1996
         
    close all
    
    % Whether to display all the params
    paramsVerbosity = 0;
    
    % Set to true if we need to compute the cone mosaic (i.e., if it does not exist)
    recomputeConeMosaic = false;
    
    % Total number of response instances to compute, say 10,000
    totalResponseInstances = 1500;
    
    % We compute instances in blocks so to enable parallel computation with
    % large cone mosaics within the limits of a machine's RAM size.
    
    % The size of the block of instances is critical for the parfor loop and 
    % must be optimized together with the number of workers. The larger the 
    % number of workers, the smaller it must be. Its value also depends on the 
    % size of the mosaic, the response length and the available RAM. 
    % For the PW96 with a 3x3 deg hex mosaic with spatially-varying 
    % cone density, and for 5 parallel workers, it cannot be more than 80 
    % for parallel computation to be realizable on a computer with 64 GB RAM.
    
    % How many response instances to be computed within each block, e.g. 80
    % with 5 workers on a 64 GB RAM system
    instancesBlockSize = 50;
        
    % Required blocks in order to generate the desired total # of response instances
    instancesBlocksNum = ceil(totalResponseInstances / instancesBlockSize);

    % Number of parallel pool workers to be spawned (>1 for parallel jobs)
    parpoolWorkersNum = 5;
    
    % How much data to save for classification
    % temporalResponseToIncludeInUnitsOfIntegrationTime = 0;     % Only peak response
    temporalResponseToIncludeInUnitsOfIntegrationTime = 17.0;  
   
    % Generate simulation params
    [sessionParams, spatialParams, temporalParams, ...
        LMSsamplingParams, chromaticDirectionParams, backgroundParams, ...
        oiParams, mosaicParams, responseSubSamplingParams] = ...
        runParamsGenerate(instancesBlocksNum, instancesBlockSize, temporalResponseToIncludeInUnitsOfIntegrationTime);
    
    % Set up the rw object for this program
    rwObject = IBIOColorDetectReadWriteBasic;
    theProgram = mfilename;
    
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
    
    % Empty the joblist
    jobDataStruct = {};
    
    % Load the jobDataStruct with all the examined conditions
    for chromaticDirectionIndex = 1:numel(chromaticDirectionParams) 
        % Load params for this chromaticDirection
        chromaticParams =  struct(...
             'backgroundxyY', backgroundParams.backgroundxyY, ...
             'coneContrasts', chromaticDirectionParams{chromaticDirectionIndex}.coneContrastUnitVector);
         
        stimulusStrengthAxis = LMSsamplingParams{chromaticDirectionIndex}.stimulusStrengthAxis;
        instancesBlocksNum = LMSsamplingParams{chromaticDirectionIndex}.instancesBlocksNum; 
        instancesBlockSize = LMSsamplingParams{chromaticDirectionIndex}.instancesBlockSize;
           
        fprintf('Chromatic direction #%d/%d. Unit vector: <%0.4f %0.4f %0.4f> \n', ...
                chromaticDirectionIndex, numel(chromaticDirectionParams), ...
                chromaticDirectionParams{chromaticDirectionIndex}.coneContrastUnitVector(1), ...
                chromaticDirectionParams{chromaticDirectionIndex}.coneContrastUnitVector(2), ...
                chromaticDirectionParams{chromaticDirectionIndex}.coneContrastUnitVector(3));
        
        for stimStrengthIndex = 1:numel(stimulusStrengthAxis)
            % Update the current stimulus strength
            chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrength = stimulusStrengthAxis(stimStrengthIndex);
            
            % Compute modulated stimulus scene
            [modulatedScene, actualStimulusStrength] = generateGaborDisplayScene(...
                spatialParams, chromaticParams, chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrength);
            
            fprintf('\tStimulus strength %d/%d, RMS cone contrast: specified = %0.5f, measured: %0.5f\n', ...
                stimStrengthIndex, numel(stimulusStrengthAxis), ...
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
                'instanceBlockIndex', 0, ...
                'instancesBlocksNum', instancesBlocksNum, ...
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
            paramsList = {mosaicParams, oiParams, spatialParams, temporalParams, backgroundParams, ...
                          LMSsamplingParams{chromaticDirectionIndex}, chromaticDirectionParams{chromaticDirectionIndex}, responseSubSamplingParams};
                  
            % Compute the oiSequence
            theOIsequence = oiSequence(oiBackground, oiModulated, ancillaryData.temporalParams.stimTimeAxis, ...
                                   ancillaryData.temporalParams.stimulusModulationFunction, 'composition', 'blend');
            %dataStruct.theOIsequence.visualize('format', 'montage');
            
            fprintf('\t\t Adding %d response blocks, each with %d response instances\n', instancesBlocksNum, instancesBlockSize);
            for instanceBlockIndex = 1:instancesBlocksNum
                
                responseData.instanceBlockIndex = instanceBlockIndex;
                
                dataStruct = struct();
                dataStruct.chromaticDirectionIndex = chromaticDirectionIndex;
                dataStruct.stimStrengthIndex = stimStrengthIndex;
                dataStruct.randomSeed = randi(1000000)+1;           % random number seed for each job      
                dataStruct.theEMpaths = [];
                dataStruct.theOIsequence =  theOIsequence;
                dataStruct.paramsList = paramsList;
                dataStruct.stimLabel = stimLabel;
                dataStruct.responseData = responseData;
                dataStruct.ancillaryData = ancillaryData;
                dataStruct.rwObject = rwObject;

                % add entry to joblist
                jobDataStruct{numel(jobDataStruct)+1} = dataStruct;
            end % for instanceBlockIndex
        end % stimStrengthIndex
    end % chromaticDirectionIndex
    
    % Each worker has its own cone mosaic
    for workerID = 1: parpoolWorkersNum
       theWorkerConeMosaic{workerID} = theConeMosaic.copy();
    end
    
    % Total number of jobs
    jobsNum = numel(jobDataStruct); 

    if (parpoolWorkersNum > 1)
        % Delete an existing parpool if one exists and start fresh
        delete(gcp('nocreate'))

        % Start parpool
        fprintf('\nOpening a local parallel pool with %d workers for dispatching %d compute jobs.\n', parpoolWorkersNum, jobsNum);
        parpool('local',parpoolWorkersNum);
    end
    
    % Loop over all jobs 
    parfor jobIndex = 1:jobsNum
        
        if (parpoolWorkersNum > 1)
            % Get the parallel pool worker ID
            t = getCurrentTask();
            workerID = t.ID;
        else
            workerID = 1;
        end
        
        % Get the data for the current condition
        theData = jobDataStruct{jobIndex};

        % Print all params
        if (paramsVerbosity > 0)
            disp(UnitTest.displayNicelyFormattedStruct(theData.ancillaryData, 'ancillaryData', '', 60));
        end
            
        % Save the ancillary data that we will need for use by the classifier preprocessing subroutine
        theData.rwObject.write('ancillaryData', theData.ancillaryData, theData.paramsList, theProgram, 'type', 'mat');
        
        % Initialize random number generator for this condition
        rng(theData.randomSeed);
        
        % Compute the emPaths for this condition
        eyeMovementsNum = computeEyeMovementsNum(theWorkerConeMosaic{workerID}.integrationTime, theData.theOIsequence);
        theData.theEMpaths = zeros(instancesBlockSize, eyeMovementsNum, 2);     
        for instanceIndex = 1:instancesBlockSize
            theData.theEMpaths(instanceIndex, :,:) = theWorkerConeMosaic{workerID}.emGenSequence(eyeMovementsNum);
        end
        if (theData.ancillaryData.mosaicParams.eyesDoNotMove)
            theData.theEMpaths= 0*theData.theEMpaths;
        end
           
        % Retrieve the responseData to be populated
        responseData = theData.responseData;
        
        % Describe what the worker is doing
        workerDescription = sprintf('Computing %d instances for job %03d/%03d [LMS=%0.3f,%0.3f,%0.3f, R=%0.4f, Block=%03d/%03d]', ...
            instancesBlockSize, jobIndex, jobsNum, ...
            responseData.stimConeContrastUnitVector(1), ...
            responseData.stimConeContrastUnitVector(2), ...
            responseData.stimConeContrastUnitVector(3), ...
            responseData.stimStrength, ...
            responseData.instanceBlockIndex, ...
            responseData.instancesBlocksNum ...
            );
        
        [responseData.absorptionsCountSequence, ...
         responseData.photoCurrentSignals] = ...
                    theWorkerConeMosaic{workerID}.computeForOISequence(theData.theOIsequence, ...
                    'emPaths', theData.theEMpaths, ...
                    'currentFlag', true, ...
                    'newNoise', true, ...
                    'workerID', workerID, ... 
                    'workDescription',  workerDescription ...
                    );
        
         responseData.absorptionsTimeAxis = theWorkerConeMosaic{workerID}.timeAxis + theData.theOIsequence.timeAxis(1);  
         responseData.photoCurrentTimeAxis = responseData.absorptionsTimeAxis;
    
        % Determine portion of the absorptions signal to be kept
        absorptionsTimeIndicesToKeep = determineTimeIndicesToKeep(responseData.absorptionsTimeAxis, ...
                theWorkerConeMosaic{workerID}.integrationTime, theData.ancillaryData.temporalParams.rampPeak, ...
                theData.ancillaryData.responseSubSamplingParams);
        
        % Determine portion of the photocurrents signal to be kept
        photocurrentsTimeIndicesToKeep = determineTimeIndicesToKeep(responseData.photoCurrentTimeAxis, ...
                theWorkerConeMosaic{workerID}.integrationTime, theData.ancillaryData.temporalParams.rampPeak, ...
                theData.ancillaryData.responseSubSamplingParams);
        
        % Remove unwanted portions of the time axes
        responseData.absorptionsTimeAxis = responseData.absorptionsTimeAxis(absorptionsTimeIndicesToKeep);
        responseData.photoCurrentTimeAxis = responseData.photoCurrentTimeAxis(photocurrentsTimeIndicesToKeep);
        
        % Remove unwanted portions of the responses
        if (isa(theWorkerConeMosaic{workerID}, 'coneMosaicHex'))
            responseData.absorptionsCountSequence = responseData.absorptionsCountSequence(:,:,absorptionsTimeIndicesToKeep);
            responseData.photoCurrentSignals = responseData.photoCurrentSignals(:,:,photocurrentsTimeIndicesToKeep);
        else
            responseData.absorptionsCountSequence = responseData.absorptionsCountSequence(:,:,:,absorptionsTimeIndicesToKeep);
            responseData.photoCurrentSignals = responseData.photoCurrentSignals(:,:,:,photocurrentsTimeIndicesToKeep);
        end

        % add the coneMosaic
        responseData.coneMosaic = theConeMosaic;
        
        % Save the computed responseData
        responseDataFile = sprintf('coneResponsesBlock%dof%d', responseData.instanceBlockIndex,instancesBlocksNum);
        rwObject.write(responseDataFile, responseData, theData.paramsList, theProgram, 'type', 'mat');
    end % all conditions parfor
    fprintf('\nAll done !\n');
end


% ------- Helper functions ------

function [sessionParams, spatialParams, temporalParams, ...
        LMSsamplingParams, chromaticDirectionParams, backgroundParams, ...
        oiParams, mosaicParams, responseSubSamplingParams] = runParamsGenerate(instancesBlocksNum, instancesBlockSize, temporalResponseToIncludeInUnitsOfIntegrationTime)
    
    sessionParams = struct(...
                 'type', 'Session', ...
          'sessionType', 'LMSplane' ...
    );

    % In P&W 1996, in the constant cycle condition, this was 10 deg (Section 2.2, p 517)
    PW96_fovDegs = 3.0;
    
    spatialParams = struct(...
                  'type', 'Spatial_v2', ...
            'spatialType', 'Gabor', ...
        'fieldOfViewDegs', PW96_fovDegs, ... 
        'cyclesPerDegree', 4.0,...
      'spatialPhaseDegs', 0, ...
       'orientationDegs', 90, ...
             'windowType', 'Gaussian', ...
        'gaussianFWHMDegs', 1.9, ...
        'viewingDistance', 0.75, ...
                    'row', 128, ...
                    'col', 128 ...
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
         'conePacking', 'hex', ...
     'fieldOfViewDegs', PW96_fovDegs,...         % nan for 1L, 1M, and 1S-cone only
    'eccentricityDegs', 0, ...  
 'spatialLMSDensities', [0.62 0.31 0.07], ...
 'integrationTimeSecs', 50/1000, ...            % 50 msec integration time
         'photonNoise', true, ...               % add Poisson noise
      'osTimeStepSecs', 10/1000, ...            % 10 milliseconds
             'osNoise', true, ...               % outer-segment noise
       'eyesDoNotMove', false ....              % normal eye movements
       );
   
  % mosaicParams.conePacking = 'rect';
  %mosaicParams.fieldOfViewDegs = PW96_fovDegs*[0.1 0.1];
   
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
    LMSsamplingParams = LMSsamplingParamsGenerate(instancesBlocksNum, instancesBlockSize);
 
    % Generate chromaticDirectionParams
    chromaticDirectionParams = chromaticDirectionParamsGenerate(LMSsamplingParams);
end
    
function LMSsamplingParams = LMSsamplingParamsGenerate(instancesBlocksNum, instancesBlockSize)
    % Add sampling params for each of the LMS directions we want to explore
    LMSsamplingParams = {};
    
    if (1==1)
        % 9 strengths along the 45,45 direction
        LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 45, ...                         % (x/y plane) (L/M modulation)
          'elevationAngle', 45, ...                        % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.1, 0.9, 9), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
      'instancesBlocksNum', instancesBlocksNum, ...
      'instancesBlockSize', instancesBlockSize ...
        );  
    
        % the null stimulus (zero stimulus strength)
        LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 0, ...                      % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                    % z-axis (S-modulation)
    'stimulusStrengthAxis', 0, ...      
      'instancesBlocksNum', instancesBlocksNum, ...
      'instancesBlockSize', instancesBlockSize ...
        );
    
    else
        
    % the null stimulus (zero stimulus strength)
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 0, ...                      % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                    % z-axis (S-modulation)
    'stimulusStrengthAxis', 0, ...      
      'instancesBlocksNum', instancesBlocksNum, ...
      'instancesBlockSize', instancesBlockSize ...
        );
    
    % L direction, specifically: cL = 1, cM = 0, cS = 0.0
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 0.0, ...                        % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                        % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.01, 0.1, 10), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
      'instancesBlocksNum', instancesBlocksNum, ...
      'instancesBlockSize', instancesBlockSize ...
        );
    
    % L+M direction, specifically: cL = 0.7071, cM = 0.7071, cS = 0.0
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 45.0, ...                        % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                        % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.01, 0.1, 10), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
      'instancesBlocksNum', instancesBlocksNum, ...
      'instancesBlockSize', instancesBlockSize ...
        );
    
    % M direction, specifically: cL = 0, cM = 1.0, cS = 0.0
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 90.0, ...                        % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                         % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.01, 0.1, 10), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
      'instancesBlocksNum', instancesBlocksNum, ...
      'instancesBlockSize', instancesBlockSize ...
        );
    
    
    % -L+M direction, specifically: cL = -0.7071, cM = 0.7071, cS = 0.0
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 135.0, ...                      % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                        % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.01, 0.1, 10), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
      'instancesBlocksNum', instancesBlocksNum, ...
      'instancesBlockSize', instancesBlockSize ...
        );  
    
    % S direction, specifically: cL = -0.0, cM = 0.0, cS = 1.0
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 0.0, ...                         % (x/y plane) (L/M modulation)
          'elevationAngle', 90.0, ...                        % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.01, 0.7, 10), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
      'instancesBlocksNum', instancesBlocksNum, ...
      'instancesBlockSize', instancesBlockSize ...
        );  
    end
    
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

    if (strcmp(mosaicParams.conePacking, 'hex')) && ~any(isnan(mosaicParams.fieldOfViewDegs))
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
        %fprintf('Optical image sequence contains %2.0f eye movements (%2.2f eye movements/oi)\n', eyeMovementsNum, eyeMovementsNumPerOpticalImage);
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

function indicesToKeep = determineTimeIndicesToKeep(timeAxis, integrationTime, rampPeak, responseSubSamplingParams)
    t = timeAxis + integrationTime/2 - rampPeak - responseSubSamplingParams.secondsToIncludeOffset;
    indicesToKeep  = find(abs(t) <= responseSubSamplingParams.secondsToInclude/2);
    if (isempty(indicesToKeep ))
       [~,indicesToKeep ] = min(abs(t-responseSubSamplingParams.secondsToInclude/2));
    end
end
             
