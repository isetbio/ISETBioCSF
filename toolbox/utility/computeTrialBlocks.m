function nParforTrials = computeTrialBlocks(ramPercentageEmployed, nTrials, coneMosaicPatternSize, coneMosaicActivePatternSize, ...
    temporalParams, spatialParams, oiRowsCols, colorModulationParams, integrationTime, displayTrialBlockPartitionDiagnostics, employStandardHostComputerResources)      
    % Determine system resources
    [numberOfWorkers, ramSizeGBytes, sizeOfDoubleInBytes] = determineSystemResources(employStandardHostComputerResources);
   
    if (numberOfWorkers == 1)
        nParforTrials = nTrials;
        return;
    end
    
    % Ensure ramPercentageEmployed is in [0.05 1]
    ramPercentageEmployed = max([0.05 ramPercentageEmployed]);
    ramSizeGBytes
    ramSizeGBytes = ramPercentageEmployed * ramSizeGBytes
    
    % Subtract RAM used by the OS
    ramUsedByOSGBytes = 1.2;
    ramSizeGBytesAvailable = ramSizeGBytes - ramUsedByOSGBytes;
    
    % Compute sizes of the large players
    if (isnan(temporalParams.windowTauInSeconds))
        [stimulusTimeAxis, ~, ~] = squareTemporalWindowCreate(temporalParams);
    else
        [stimulusTimeAxis, ~, ~] = gaussianTemporalWindowCreate(temporalParams);
    end
    
    eyeMovementsNumPerOpticalImage = max([1 ceil(temporalParams.stimulusSamplingIntervalInSeconds/integrationTime)]);
    emPathLength = eyeMovementsNumPerOpticalImage*numel(stimulusTimeAxis);
    
    wavelengths = floor((colorModulationParams.endWl-colorModulationParams.startWl)/colorModulationParams.deltaWl);
    opticalImageSize = 2.5*oiRowsCols(1)*oiRowsCols(2)*wavelengths;

    % estimate sizes of the various matrices used
    trialBlockSize = floor(nTrials/numberOfWorkers);
    totalMemoryPerWorker = computeTotalMemoryPerWorker();
    totalMemoryUsed = totalMemoryPerWorker * numberOfWorkers;
    
    allowedRAMcompression = 1.0;
    while (totalMemoryUsed > allowedRAMcompression*ramSizeGBytesAvailable)
        [totalMemoryUsed  allowedRAMcompression*ramSizeGBytesAvailable trialBlockSize]
        trialBlockSize = trialBlockSize-1;
        totalMemoryPerWorker = computeTotalMemoryPerWorker();
        totalMemoryUsed = totalMemoryPerWorker * numberOfWorkers;
    end
    fprintf('trialBlockSize: %d, total memory used: %f\n', trialBlockSize, totalMemoryUsed);
    
    trialBlockSize = min([nTrials max([1 trialBlockSize])]);
    nParforTrialBlocks = ceil(nTrials / trialBlockSize);
    totalMemoryPerWorker = computeTotalMemoryPerWorker();
    
    if (nParforTrialBlocks < numberOfWorkers)
        nParforTrialBlocks = numberOfWorkers;
        trialBlockSize = max([1 floor(nTrials/nParforTrialBlocks)]);
        totalMemoryPerWorker = computeTotalMemoryPerWorker();
    end
    
    totalMemoryUsed = numberOfWorkers * totalMemoryPerWorker;
    
    if (totalMemoryUsed < 0.01*ramSizeGBytesAvailable/numberOfWorkers)
         % Just use one processor - faster most of the times
         nParforTrials(1) = nTrials;
         nParforTrialBlocks = 1;
         trialBlockSize = nTrials;
    elseif (totalMemoryUsed < ramSizeGBytesAvailable/numberOfWorkers)
          nParforTrialBlocks = numberOfWorkers;
          trialBlockSize = floor(nTrials/nParforTrialBlocks);
          for kk = 1:nParforTrialBlocks-1
              nParforTrials(kk) = trialBlockSize;
          end
          nParforTrials(nParforTrialBlocks) = nTrials - trialBlockSize*nParforTrialBlocks;
          totalMemoryPerWorker = computeTotalMemoryPerWorker();
          totalMemoryUsed = totalMemoryPerWorker * numberOfWorkers;
    end

     
    % Last block of trials
    nParforTrialsTmp = ones(1, nParforTrialBlocks) * trialBlockSize;
    remainingTrials = nTrials - sum(nParforTrialsTmp);
    if (remainingTrials > 0) 
        if remainingTrials > round(trialBlockSize*0.2)
            nParforTrialsTmp(end+1) = remainingTrials;
        else
            nParforTrialsTmp(end) = nParforTrialsTmp(end) + remainingTrials;
        end
    end
    
    % Check that # of trials is correct
    accumTrials = 0;
    k = 1; keepGoing = true;
    while (k <= numel(nParforTrialsTmp)) && (keepGoing)
        if (accumTrials+nParforTrialsTmp(k) > nTrials)
            if (nTrials > accumTrials)
                nParforTrials(k) = nTrials - accumTrials;
            end
            keepGoing = false;
        else
            nParforTrials(k) = nParforTrialsTmp(k);
            accumTrials = accumTrials+nParforTrials(k);
            k = k + 1;
        end
    end

    if (sum(nParforTrials) ~= nTrials)
        nParforTrials
        nTrials
        trialBlockSize
        nParforTrialBlocks
        error('Error in logic of trial partitioning.')
    end
    
    fprintf('<strong> %d workers, system RAM = %2.1fGBytes </strong> \n', numberOfWorkers, ramSizeGBytes), ...
    fprintf('<strong> %d trials partitioned in %d blocks, each with %d trials (last has %d trials) </strong>\n', nTrials, numel(nParforTrials), nParforTrials(1), nParforTrials(end));
    fprintf('<strong> RAM used : %2.1f GBytes (per worker), total: %2.1f GBytes </strong> \n\n', totalMemoryPerWorker, totalMemoryUsed);
    
    % Nested function 
    function totalMemoryPerWorker = computeTotalMemoryPerWorker()
        totalMemoryPerWorker = (opticalImageSize + 2*coneMosaicPatternSize*trialBlockSize*eyeMovementsNumPerOpticalImage + coneMosaicActivePatternSize*emPathLength*trialBlockSize)*sizeOfDoubleInBytes/(1024^3);
    end

end

