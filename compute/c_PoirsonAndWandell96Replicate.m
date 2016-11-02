function c_PoirsonAndWandell96Replicate
% c_PoirsonAndWandell96Replicate
%
% Compute color detection thresholds to replicate the Poirson & Wandell 1996
         
    close all
    
    % Set up the rw object for this program
    rwObject = IBIOColorDetectReadWriteBasic;
    theProgram = mfilename;

    sessionParams = struct(...
                 'type', 'Session', ...
          'sessionType', 'LMSplane' ...
    );

    spatialParams = struct(...
                  'type', 'Spatial_v2', ...
            'spatialType', 'Gabor', ...
        'fieldOfViewDegs', 5.0, ...      In P&W 1996, in the constant cycle condition, this was 10 deg (Section 2.2, p 517)
        'cyclesPerDegree', 2.0,...
      'spatialPhaseInDeg', 0, ...
       'orientationInDeg', 0, ...
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
             'rampDurationSecs', 50.0 * 1.0/CRTrefreshRate, ... 
                  'rampTauSecs', 165/1000, ...
                     'rampPeak', 0/1000 ...
        );
    
    
    % Define the stimulus temporal modulation function
    temporalParams.stimTimeAxis = 0:temporalParams.stimulusSamplingIntervalSecs:temporalParams.rampDurationSecs;
    temporalParams.stimTimeAxis =  temporalParams.stimTimeAxis - mean(temporalParams.stimTimeAxis);
    % Compute modulation function
    temporalParams.stimulusModulationFunction = exp(-0.5*((temporalParams.stimTimeAxis-temporalParams.rampPeak)/temporalParams.rampTauSecs).^2);
    
    % Optics params
    oiParams = oiParamsGenerate();
    oiParams.type = 'Optics_v2';
    oiParams.fieldOfViewDegs = spatialParams.fieldOfViewDegs*1.2;
    
    mosaicParams = struct(...
                'type', 'Mosaic_v2', ...
          'mosaicType', 'RECT', ...
     'fieldOfViewDegs', 0.1*spatialParams.fieldOfViewDegs,...         % nan for 1L, 1M, and 1S-cone only
    'eccentricityDegs', 0, ...  
 'spatialLMSDensities', [0.6 0.3 0.1], ...
 'integrationTimeSecs', 50/1000, ...        % 50 msec integration time
         'photonNoise', true, ...           % add Poisson noise
      'osTimeStepSecs', 5/1000, ...         % 5 milliseconds
             'osNoise', true, ...           % outer-segment noise
       'eyesDoNotMove', false ....          % normal eye movements
       );
    
   responseSubSamplingParams = struct(...
                     'type', 'ResponseSubsampling', ...
          'secondsToInclude', mosaicParams.integrationTimeSecs*0.0, ...    % temporal subsampling: only keep responses within this time period around 'temporalParams.rampPeak'
    'secondsToIncludeOffset', 0 ...                                        % and this offset
        );
   
    % In the constant cycle condition, the background was xyY= 0.38, 0.39, 536.2 cd/m2
    % Also they say they placed a uniform field to the screen to increase the contrast resolutions (page 517, Experimental Aparratus section)
    xyY =  [0.38 0.39 536.2];
    backgroundChromaticParams = struct(...
                     'type', 'Background_v2', ...
                   'device', 'Monitor', ...
            'backgroundxyY', xyY, ...
            'coneContrasts', [0.0 0.0 0.0] ...
        );
    
    % Add sampling params for each of the LMS directions we want to explore
    LMSsamplingParams = {};
    
    % the null stimulus (zero stimulus strength)
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 0, ...                      % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                    % z-axis (S-modulation)
    'stimulusStrengthAxis', 0, ...      
            'instancesNum', 100 ...
        );
    
    % -L+M direction, specifically: cL = -0.7071, cM = 0.7071, cS = 0.0
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                   'type', 'LMSsampling', ...
            'azimuthAngle', 135.0, ...                      % (x/y plane) (L/M modulation)
          'elevationAngle', 0.0, ...                        % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.01, 0.1, 5), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
            'instancesNum', 100 ...
        );
    
    % L+M+S, specifically: cL = 0.5, cM = 0.5, cS = 0.7071
    LMSsamplingParams{numel(LMSsamplingParams)+1} = struct(...
                    'type', 'LMSsampling', ...
            'azimuthAngle', 45.0, ...                       % (x/y plane) (L/M modulation)
          'elevationAngle', 45.0, ...                       % z-axis (S-modulation)
    'stimulusStrengthAxis', linspace(0.01, 0.1, 5), ...     % linspace(min, max, nLevels) or logspace(log10(min), log10(max), nLevels)
            'instancesNum', 100 ...
        );
    
    
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
    
    % Generate the background scene
    [backgroundScene, ~] = generateGaborDisplayScene(spatialParams, backgroundChromaticParams, 0.0);
    
    % Generate custom optics
    theOI = colorDetectOpticalImageConstruct(oiParams);
    
    % Compute the background OI
    oiBackground = theOI;
    oiBackground = oiCompute(oiBackground, backgroundScene);
  
    % Loop over the examined chromatic directions
    for chromaticDirectionIndex = 1:numel(chromaticDirectionParams)  
        % Load params for this chromaticDirection
        stimulusChromaticParams =  struct(...
             'backgroundxyY', backgroundChromaticParams.backgroundxyY, ...
             'coneContrasts', chromaticDirectionParams{chromaticDirectionIndex}.coneContrastUnitVector);
         
         stimulusStrengthAxis = LMSsamplingParams{chromaticDirectionIndex}.stimulusStrengthAxis;
         instancesNum = LMSsamplingParams{chromaticDirectionIndex}.instancesNum;
                  
        % Compute responses for all stimulus strengths
        for stimStrengthIndex = 1:numel(stimulusStrengthAxis)
            % Update the current stimulus strength
            chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrength = stimulusStrengthAxis(stimStrengthIndex);
            
            % Compute modulated stimulus scene
            [modulatedScene, actualStimulusStrength] = generateGaborDisplayScene(spatialParams, stimulusChromaticParams, chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrength);
            fprintf('Stimulus %d/%d - Strength (RMS cone contrast): specified = %2.3f, measured: %2.3f\n', stimStrengthIndex, numel(stimulusStrengthAxis), stimulusStrengthAxis(stimStrengthIndex), actualStimulusStrength);
            
            % Compute modulated OI
            oiModulated = theOI;
            oiModulated = oiCompute(oiModulated, modulatedScene);
    
            % Generate the sequence of optical images representing the ramping of the stimulus
            theOIsequence = oiSequence(oiBackground, oiModulated, temporalParams.stimTimeAxis, temporalParams.stimulusModulationFunction, 'composition', 'blend');
            theOIsequence.visualize();
    
            if ((chromaticDirectionIndex == 1) && (stimStrengthIndex == 1))
                % Generate the cone mosaic
                [theConeMosaic, eyeMovementsNum] = coneMosaicGenerate(mosaicParams, theOIsequence); 
            end
            
            % Generate stimulus label
            stimLabel = sprintf('AzimuthAngle=%2.1fDeg, ElevationAngle=%2.1fDeg, StimStrength=%2.3f', ...
                LMSsamplingParams{chromaticDirectionIndex}.azimuthAngle, ...
                LMSsamplingParams{chromaticDirectionIndex}.elevationAngle, ...
                chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrength);
            % Generate response storage struct
            stimData = struct(...
                'stimLabel', stimLabel, ...
                'stimConeContrastUnitVector', chromaticDirectionParams{chromaticDirectionIndex}.coneContrastUnitVector, ...
                'stimStrength',  chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrength, ...
                'absorptionsCountSequence', [], ...
                'absorptionsTimeAxis', [], ... 
                'photoCurrentSignals', [], ... 
                'photoCurrentTimeAxis', []);
            
            ancillaryData = struct(...
                'sessionParams', sessionParams, ...
                'oiParams', oiParams, ...
                'mosaicParams', mosaicParams, ...
                'spatialParams', spatialParams, ...
                'temporalParams', temporalParams, ...
                'backgroundChromaticParams',  backgroundChromaticParams, ...
                'LMSsamplingParams',  LMSsamplingParams{chromaticDirectionIndex}, ...
                'chromaticDirectionParams', chromaticDirectionParams{chromaticDirectionIndex}, ...
                'responseSubSamplingParams', responseSubSamplingParams ...
            );
            
            % Print all params
            disp(UnitTest.displayNicelyFormattedStruct(ancillaryData, 'ancillaryData', '', 60));
            
            fprintf('Computing emPaths\n');
            % Compute the emPaths for all response instances
            for instanceIndex = 1:instancesNum
                emPaths(instanceIndex, :,:) = theConeMosaic.emGenSequence(eyeMovementsNum);
            end
            if (mosaicParams.eyesDoNotMove)
                emPaths = 0*emPaths;
            end
            
            fprintf('%s\nComputing %d response instances for a %d x %d cone mosaic with %d eye movements scanning a %d-frame stimulus sequence.\n', ...
                stimLabel, instancesNum, theConeMosaic.mosaicSize(1), theConeMosaic.mosaicSize(2), eyeMovementsNum, theOIsequence.length);
            % Compute absorptions and photocurrents for all response instances
            tic
            [stimData.absorptionsCountSequence, stimData.absorptionsTimeAxis, stimData.photoCurrentSignals, stimData.photoCurrentTimeAxis] = ...
                    theConeMosaic.computeForOISequence(theOIsequence, ...
                    'emPaths', emPaths, ...
                    'currentFlag', true, ...
                    'newNoise', true ...
                    );
            fprintf('Response computation took %2.2f minutes\n', toc/60);
    
            % Temporal response subsampling for absorptions
            t = stimData.absorptionsTimeAxis + theConeMosaic.integrationTime/2 - temporalParams.rampPeak - responseSubSamplingParams.secondsToIncludeOffset;
            timeIndicesToKeep = find(abs(t) <= responseSubSamplingParams.secondsToInclude/2);
            if (isempty(timeIndicesToKeep))
                [~,timeIndicesToKeep] = min(abs(t-responseSubSamplingParams.secondsToInclude/2));
            end
            
            % remove unwanted portions of the absorption responses
            stimData.absorptionsTimeAxis = stimData.absorptionsTimeAxis(timeIndicesToKeep);
            stimData.absorptionsCountSequence = stimData.absorptionsCountSequence(:,:,:,timeIndicesToKeep);
            
            % Temporal response subsampling for photocurrents
            t = stimData.photoCurrentTimeAxis + theConeMosaic.integrationTime/2 - temporalParams.rampPeak -responseSubSamplingParams.secondsToIncludeOffset;
            timeIndicesToKeep = find(abs(t) <= responseSubSamplingParams.secondsToInclude/2);
            if (isempty(timeIndicesToKeep))
                [~,timeIndicesToKeep] = min(abs(t-responseSubSamplingParams.secondsToInclude/2));
            end
            
            % remove unwanted portions of the photocurrent responses
            stimData.photoCurrentTimeAxis = stimData.photoCurrentTimeAxis(timeIndicesToKeep);
            stimData.photoCurrentSignals = stimData.photoCurrentSignals(:,:,:,timeIndicesToKeep);
            
            %
            % Save responses and parameters for this stimStrengthIndex and chromaticDirectionIndex
            paramsList = {sessionParams, oiParams, mosaicParams, spatialParams, temporalParams, backgroundChromaticParams, ...
                          LMSsamplingParams{chromaticDirectionIndex}, chromaticDirectionParams{chromaticDirectionIndex}, responseSubSamplingParams};
            rwObject.write('coneResponses', stimData, paramsList, theProgram, 'type', 'mat');
            
            % Save other data we need for use by the classifier preprocessing subroutine
            rwObject.write('ancillaryData', ancillaryData, paramsList, theProgram, 'type', 'mat');
    
            % Visualize oisequence with eye movement sequence
            % theOIsequence.visualizeWithEyeMovementSequence(absorptionsTimeAxis);
    
            if ((chromaticDirectionIndex == 1) && (stimStrengthIndex == 1))
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
            end
            
          % Plot all response instances for the center-most L,M, and S cone
          instancesToPlot = 1:instancesNum;
          coneStride = 1e12;
          plotAbsorptionsCountTimeSeries(...
              stimData.absorptionsTimeAxis, ...
              stimData.absorptionsCountSequence(instancesToPlot,:,:,:), ...
              theConeMosaic.integrationTime, ...
              temporalParams.stimTimeAxis, temporalParams.stimulusModulationFunction, ...
              iL(idxL), iM(idxM), iS(idxS), coneStride, chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrength);
         
          % Plot some L,M, and S cone response for the first intance 
          instancesToPlot = 1;
          coneStride = 20;  % how many cones to skip over
          plotAbsorptionsCountTimeSeries(...
              stimData.absorptionsTimeAxis, ...
              stimData.absorptionsCountSequence(instancesToPlot,:,:,:), ...
              theConeMosaic.integrationTime, ...
              temporalParams.stimTimeAxis, temporalParams.stimulusModulationFunction, ...
              iL, iM, iS,  coneStride, chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrength);
          clear 'tmp'
        end % for contrastIndex
    end %  for chromaticDirectionIndex
    
end

function plotAbsorptionsCountTimeSeries(timeAxis, absorptions, integrationTime, stimTimeAxis, stimulusModulationFunction, iL, iM, iS,  coneStride, stimulusStrength)
    
    instancesPlotted = size(absorptions,1);
    timePointsPlotted = numel(timeAxis);
    
    minAbsorptions = min(absorptions(:));
    maxAbsorptions = max(absorptions(:));
    
    % Extract cone responses by type
    iL = iL(1:coneStride:end);
    iM = iM(1:coneStride:end);
    iS = iS(1:coneStride:end);
    if (timePointsPlotted == 1)
        absorptionsByConeType{1} = reshape(absorptions(:,iL), [instancesPlotted*numel(iL) 1]);
        absorptionsByConeType{2} = reshape(absorptions(:,iM), [instancesPlotted*numel(iM) 1]);
        absorptionsByConeType{3} = reshape(absorptions(:,iS), [instancesPlotted*numel(iS) 1]);
    else
        absorptions = permute(reshape(absorptions, [instancesPlotted size(absorptions,2)*size(absorptions,3) timePointsPlotted]), [1 3 2]);
        absorptionsByConeType{1} = reshape(permute(absorptions(:,:,iL), [1 3 2]), [instancesPlotted*numel(iL) timePointsPlotted]);
        absorptionsByConeType{2} = reshape(permute(absorptions(:,:,iM), [1 3 2]), [instancesPlotted*numel(iM) timePointsPlotted]);
        absorptionsByConeType{3} = reshape(permute(absorptions(:,:,iS), [1 3 2]), [instancesPlotted*numel(iS) size(absorptions,2)]);
    end
    colors = [1 0 0; 0 1.0 0; 0 0.8 1];
    plotBackgroundColor = [0.1 0.1 0.1];
    barOpacity = 0.1;
    
    hFig = figure(); clf;
    set(hFig, 'Color', [0 0 0], 'Position', [10 10 650 900]);
    subplot('Position', [0.08 0.08 0.85 0.87]);
    
    yyaxis left
    hold on
    % Identify stimulus presentation times
    for k = 1:numel(stimTimeAxis)
        plot(stimTimeAxis(k)*[1 1]*1000, [minAbsorptions maxAbsorptions], 'k-', 'Color', [0.5 0.5 0.5]);
    end
        
    dt = integrationTime;
    % Plot absorptions
    for coneType = 1:3
        for tIndex = 1:numel(timeAxis)
            quantaAtThisTimeBin = squeeze(absorptionsByConeType{coneType}(:,tIndex));
            plot([timeAxis(tIndex) timeAxis(tIndex)+dt]*1000, [quantaAtThisTimeBin(:) quantaAtThisTimeBin(:)], '-', 'LineWidth', 1.5, 'Color', [colors(coneType,:) barOpacity]);
        end
    end
    
    
    box on;
    xlabel('time (ms)', 'FontSize', 14);
    ylabel(sprintf('# of absorptions per %2.1f ms',dt*1000), 'FontSize', 14);
    set(gca, 'XLim', [min([stimTimeAxis(1) timeAxis(1)-min([integrationTime stimTimeAxis(2)-stimTimeAxis(1)])]) max([stimTimeAxis(end) timeAxis(end)+dt])]*1000, 'YLim', [minAbsorptions maxAbsorptions], ...
             'Color', plotBackgroundColor, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], ...
             'FontSize', 12);
         
    % Plot the stimulus modulation on the right time axis
    yyaxis right
    dt = stimTimeAxis(2)-stimTimeAxis(1);
    for tIndex = 1:numel(stimTimeAxis)
        plot([stimTimeAxis(tIndex) stimTimeAxis(tIndex)+dt]*1000, stimulusModulationFunction(tIndex)*[1 1], '-', 'LineWidth', 1.5, 'Color', 'y');
        if (tIndex < numel(stimTimeAxis))
            plot(stimTimeAxis(tIndex+1)*[1 1]*1000, [stimulusModulationFunction(tIndex) stimulusModulationFunction(tIndex+1)], '-', 'LineWidth', 1.5, 'Color', 'y');
        end
    end
    set(gca, 'YColor', 'y', 'FontSize', 12);
    ylabel('stimulus modulation', 'FontSize', 14);
     
    title(sprintf('[Stim. strength: %2.3f] Instances plotted: %d, Cones responses plotted: %d(L), %d(M), %d(S)', stimulusStrength, instancesPlotted, numel(iL), numel(iM), numel(iS)), 'Color', [1 1 1], 'FontSize',14);
    drawnow;
end


function [theConeMosaic, eyeMovementsNum] = coneMosaicGenerate(mosaicParams, theOIsequence)
    % Default human mosaic
    theConeMosaic = coneMosaic;
    
    % Adjust size
    if isnan(mosaicParams.fieldOfViewDegs)
        % Generate a human cone mosaic with 1L, 1M and 1S cone
        theConeMosaic.rows = 1;
        theConeMosaic.cols = 3;
        theConeMosaic.pattern = [2 3 4];
    else
        theConeMosaic.setSizeToFOV(mosaicParams.fieldOfViewDegs);
    end
    
    % Set the LMS spatial densities
    theConeMosaic.spatialDensity = [0 mosaicParams.spatialLMSDensities]';
    
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

    % Generate eye movement sequence for all oi's
    stimulusSamplingInterval = theOIsequence.oiTimeAxis(2)-theOIsequence.oiTimeAxis(1);
    eyeMovementsNumPerOpticalImage = stimulusSamplingInterval/theConeMosaic.integrationTime;
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