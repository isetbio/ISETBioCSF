function testPhotocurrentModel

    recomputeResponses = true;
    
    % Pulse duration in seconds
    
    %vParams.pulseDurationSeconds = 20/1000;
    vParams.pulseDurationSeconds = 100/1000;
    %vParams.pulseDurationSeconds = 200/1000;
    %vParams.pulseDurationSeconds = 400/1000;
    
    % Data filename
    dataFileName = sprintf('results_%dmsec.mat', vParams.pulseDurationSeconds*1000);
        
    % Recompute responses    
    if (recomputeResponses) 
        fprintf('Will save to %s\n', dataFileName);
        
        % Constant params
        cParams.spontaneousIsomerizationRate = 0; %  R*/c/s
        cParams.eccentricity  = 'foveal';
        cParams.useDefaultImplementation = true;
        cParams.noisyInstancesNum = 1000;
        
        % Varied Params
        % Pulse contrast
        vParams.weberContrast = -95/100;
        
        % Adaptation photon rate in R*/c/s
        vParams.adaptationPhotonRate = 12000;
        
        % The photonIntegrationTime affects the SNR of the cone
        % excitations. The longer it is the higher the SNR. 
        % Here we setting it to the pulse duration.
        vParams.photonIntegrationTime = vParams.pulseDurationSeconds;

        % Examined contrast levels
        contrastLevels = [0.02 0.04 0.08 0.16 0.32 0.64 1.0];
        contrastLevels = cat(2, -contrastLevels,  contrastLevels);

        % Examined adaptation levels
        % prctile               1%      5%     10%     25%     50%     75%     90%    95%      99% 
        %adaptationLevels = [0.0317  0.1434  0.2882  0.7143  1.4192  2.1320 2.5554  2.6961  2.8136]*1e4; % for 2 deg mosaic
        adaptationLevels = [250 500 1000 2000 4000 8000 16000]; 

        % Preallocate memory
        d = cell(1,numel(adaptationLevels));
        
        % Run all conditions
        for iAdaptationIndex = 1:numel(adaptationLevels)
            dd = cell(1,numel(contrastLevels));
            theAdaptationLevel = adaptationLevels(iAdaptationIndex);
            tic
            parfor iContrastIndex = 1:numel(contrastLevels)
                fprintf('bkgnd = %d (of %d); c = %d (of %d)\n', iAdaptationIndex, numel(adaptationLevels), iContrastIndex, numel(contrastLevels));
                theVParams = vParams;
                theVParams.weberContrast = contrastLevels(iContrastIndex);
                theVParams.adaptationPhotonRate = theAdaptationLevel;
                dd{iContrastIndex} = runSimulation(theVParams, cParams); 
            end
            fprintf('Finished in %2.1f seconds\n', toc);
            d{iAdaptationIndex} = dd;
        end

        size(d)
        save(dataFileName, 'd', 'contrastLevels', 'adaptationLevels', '-v7.3');
        fprintf('Saved data to %s\n', dataFileName);
        
    else
        % Load previously computed responses
        fprintf('Will load from %s\n', dataFileName);
        load(dataFileName, 'd', 'contrastLevels', 'adaptationLevels');       
    end
    
    for iAdaptationIndex = 1:numel(adaptationLevels)
        dd = d{iAdaptationIndex};
        for iContrastIndex = 1:numel(contrastLevels)
            theConeExcitationSNR(iAdaptationIndex, iContrastIndex) = dd{iContrastIndex}.theConeExcitationSNR;
            thePhotoCurrentSNR(iAdaptationIndex, iContrastIndex) = dd{iContrastIndex}.thePhotoCurrentSNR;
        end
    end  
        
    
    % Plot results as a function of adaptation photon rate
    % Below lims and ticks for 100 msec pulse
    SNRLims = [0.04 100];
    SNRTicks = [0.03 0.1 0.3 1 3 10 30 100 300];
    SNRratioLims = [0.01 1.0];
    SNRratioTicks = [0:0.2:1.0];
    figNo = 4;
    plotSNRanalysis(adaptationLevels, theConeExcitationSNR, thePhotoCurrentSNR, contrastLevels,  ...
        SNRLims, SNRTicks, SNRratioLims, SNRratioTicks, vParams.pulseDurationSeconds, figNo, ...
        'adaptationLevel');
    
    figNo = 400;
    plotSNRanalysis(adaptationLevels, theConeExcitationSNR, thePhotoCurrentSNR, contrastLevels,  ...
        SNRLims, SNRTicks, SNRratioLims, SNRratioTicks, vParams.pulseDurationSeconds, figNo, ...
        'contrastLevel');
    
    
    coneExcitationRange = [0 3000];  % or enter [] for automatic scaling
    
    % Plot the responses
    combinePolarities = true;
    showSNR = true;
    showSNRcomponents = ~true;
    photoCurrentRange = [-80 0];
    
    adaptationLevelsVisualized = 8000; %[2000 4000 8000 16000 24000];  % adaptationLevels
    contrastLevelsVisualized = [0.16 0.32 0.64 1.0];  % contrastLevels (specify only positive)
    
    if (combinePolarities)
        % Plot separate responses
        for iAdaptationIndex = 1:numel(adaptationLevelsVisualized)
            dd = d{iAdaptationIndex};
            theAdaptationIndex = find(adaptationLevels == adaptationLevelsVisualized(iAdaptationIndex));
            if isempty(theAdaptationIndex)
                error('Could not find data for adaptation level %2.1f R*/c/s', adaptationLevelsVisualized(iAdaptationIndex));
            end
            for iContrastIndex = 1:numel(contrastLevelsVisualized)
                theContrastIndex = find(contrastLevels == -contrastLevelsVisualized(iContrastIndex));
                if isempty(theAdaptationIndex)
                    error('Could not find data for contrast level %2.1f R*/c/s', -contrastLevelsVisualized(iContrastIndex));
                end
                dNegativePolarity = dd{theContrastIndex};
                dPositivePolarity = dd{theContrastIndex +numel(contrastLevels)/2};
                figNo = 1;
                plotResponses(dNegativePolarity.modelResponse, dPositivePolarity.modelResponse,adaptationLevels(theAdaptationIndex), ...
                    contrastLevels(theContrastIndex), vParams.pulseDurationSeconds, photoCurrentRange, showSNR, showSNRcomponents, coneExcitationRange, figNo);
            end
        end
    else
        % Plot separate responses
        for iAdaptationIndex = 1:numel(adaptationLevels)
            dd = d{iAdaptationIndex};
            for iContrastIndex = 1:numel(contrastLevels)
                dSelected = dd{iContrastIndex};
                figNo = 1;
                plotResponses(dSelected.modelResponse, [], adaptationLevels(iAdaptationIndex), ...
                    contrastLevels(iContrastIndex), vParams.pulseDurationSeconds, photoCurrentRange,  showSNR, showSNRcomponents, coneExcitationRange, figNo);
            end
        end
    end
    
    
%     plotModelResponse(dSelect.modelResponse, ...
%         dSelect.theConeExcitationSNR, ...
%         dSelect.thePhotoCurrentSNR, ...
%         dSelect.noiseEstimationLatency, ...
%         dSelect.coneExcitationModulationPeak, ...
%         dSelect.coneExcitationPhotocurrentNoiseSigma, ...
%         dSelect.photocurrentModulationPeak, ...
%         dSelect.photocurrentNoiseSigma); 
end


    
