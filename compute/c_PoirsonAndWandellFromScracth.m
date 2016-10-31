function c_PoirsonAndWandellFromScractch
% c_PoirsonAndWandellFromScratch(varargin)
%
% Compute color detection thresholds to replicate the Poirson & Wandell 1996
         
    close all
    spatialParams = struct(...
        'fieldOfViewDegs', 5.0, ...      In P&W 1996, in the constant cycle condition, this was 10 deg (Section 2.2, p 517)
        'cyclesPerDegree', 2.0,...
      'spatialPhaseInDeg', 0, ...
             'windowType', 'Gaussian', ...
        'gaussianFWHMDegs', 1.9, ...
        'viewingDistance', 0.75, ...
                    'row', 256, ...
                    'col', 256 ...
        );
    
    CRTrefreshRate = 87;
    temporalParams = struct(...
        'stimSamplingInterval', 1.0/CRTrefreshRate, ...
                 'envelopeTau', 165/1000, ...
               'envelopeWidth', 30.0/CRTrefreshRate, ... 
                'envelopePeak', 0/1000 ...
        );
    
    mosaicParams = struct(...
                'size', 0.1*spatialParams.fieldOfViewDegs,...         % nan for 1L, 1M, and 1S-cone only
     'integrationTime', 50/1000, ...        % 50 msec integration time
         'photonNoise', true, ...           % add Poisson noise
          'osTimeStep', 5/1000, ...         % 5 milliseconds
             'osNoise', true ...            % outer-segment noise
       );
    
    
    % In the constant cycle condition, the background was xyY= 0.38, 0.39, 536.2 cd/m2
    % Also they say they placed a uniform field to the screen to increase the contrast resolutionso (page 517, Experimental Aparratus section)
    xyY =  [0.38 0.39 536.2];
    backgroundChromaticParams = struct(...
            'backgroundxyY', xyY, ...
            'coneContrasts', [0.0 0.0 0.0] ...
        );
    
    % Chromatic direction angles: azimuth (LM plane) and elevation 
    % -L+M, specifically: cL = -0.5, cM = 0.5, 
    azimuthAnglesTested = [135]; 
    elevationAnglesTested = [0];
    
    % L+M+S, specifically: cL = 0.5, cM = 0.5, cS = 0.7071
    azimuthsTested = [45]; 
    elevationsTested = [45];
    
    chromaticDirectionParams = {};
    for azimuthIndex = 1:numel(azimuthAnglesTested)
        % Get angle on LM plane
        azimuthAngle = azimuthAnglesTested(azimuthIndex);
        
        for elevationsIndex = 1:numel(elevationsTested)
            % Get elevation angle from LM plane
            elevationAngle = elevationAnglesTested(elevationsIndex);
            
            chromaticDirectionParams{numel(chromaticDirectionParams)+1} = struct(...
                      'azimuthAngle', azimuthAngle, ...
                    'elevationAngle', elevationAngle, ...
                      'instancesNum', 100, ...
              'stimulusStrengthsNum', 5, ...
               'minStimulusStrength', 0.01, ...
               'maxStimulusStrength', 0.15, ...
             'stimulusStrengthScale', 'linear');
        end  % elevationIndex
    end  % azimuthIndex

    
    % Define the stimulus temporal modulation function
    stimTimeAxis = -temporalParams.envelopeWidth:temporalParams.stimSamplingInterval:temporalParams.envelopeWidth;
    % monophasic modulation function
    stimulusModulationFunction = exp(-0.5*((stimTimeAxis-temporalParams.envelopePeak)/temporalParams.envelopeTau).^2);
    
    
    % Generate the background scene
    [backgroundScene, ~] = generateGaborDisplayScene(spatialParams, backgroundChromaticParams, 0.0);
    
    % Generate custom optics
    oiParams = struct(...
        'fieldOfViewDegs', 5.0 ...
        );
    theOI = oiGenerate(oiParams);
    
    % Compute the background OI
    oiBackground = theOI;
    oiBackground = oiCompute(oiBackground, backgroundScene);

    
    % Loop over the examined chromatic directions
    for chromaticDirectionIndex = 1:numel(chromaticDirectionParams)  
        % Compute cone contrasts from azimuth and elevation
        [cL, cM, cS] = sph2cart(chromaticDirectionParams{chromaticDirectionIndex}.azimuthAngle/180*pi, chromaticDirectionParams{chromaticDirectionIndex}.elevationAngle/180*pi, 1);
        
        % Normalize to unity RMS cone contrast (stimulus strength)
        s = [cL, cM, cS];
        s = s / norm(s);
        
        stimulusChromaticParams =  struct(...
                        'backgroundxyY', backgroundChromaticParams.backgroundxyY, ...
                        'coneContrasts', s);
        
        if (strcmp(chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrengthScale, 'linear'))
                nominalStimulusStrengths = linspace(...
                    chromaticDirectionParams{chromaticDirectionIndex}.minStimulusStrength, ...
                    chromaticDirectionParams{chromaticDirectionIndex}.maxStimulusStrength, ...
                    chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrengthsNum);
        elseif (strcmp(chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrengthScale, 'log'))
                nominalStimulusStrengths  = logspace(...
                    log10(chromaticDirectionParams{chromaticDirectionIndex}.minStimulusStrength), ...
                    log10(chromaticDirectionParams{chromaticDirectionIndex}.maxStimulusStrength), ...
                    chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrengthsNum);
        else
            error('Unknown stimulus strength scale: ''%s''\n', chromaticDirectionParams{chromaticDirectionIndex}.stimulusStrengthScale);
        end
            
        instancesNum = chromaticDirectionParams{chromaticDirectionIndex}.instancesNum;
        
        % Compute responses for all stimulus strengths
        for stimStrengthIndex = 1:numel(nominalStimulusStrengths)
            
            % Compute modulated stimulus scene
            [modulatedScene, actualStimulusStrength] = generateGaborDisplayScene(spatialParams, stimulusChromaticParams, nominalStimulusStrengths(stimStrengthIndex));
            fprintf('Strength (RMS cone contrast): specified = %2.3f, measured: %2.3f\n', nominalStimulusStrengths(stimStrengthIndex), actualStimulusStrength);
            
            % Compute modulated OI
            oiModulated = theOI;
            oiModulated = oiCompute(oiModulated, modulatedScene);
    
            % Generate the sequence of optical images representing the ramping of the stimulus
            theOIsequence = oiSequence(oiBackground, oiModulated, stimTimeAxis, stimulusModulationFunction, 'composition', 'blend');
            %theOIsequence.visualize();
    
            if ((chromaticDirectionIndex == 1) && (stimStrengthIndex == 1))
                % Generate the cone mosaic
                [theConeMosaic, eyeMovementsNum] = coneMosaicGenerate(mosaicParams, theOIsequence); 
            end
            
            fprintf('Computing emPaths\n');
            % Compute the emPaths for all response instances
            for instanceIndex = 1:chromaticDirectionParams{chromaticDirectionIndex}.instancesNum
                emPaths(instanceIndex, :,:) = theConeMosaic.emGenSequence(eyeMovementsNum);
            end
    
            % Compute absorptions and photocurrents for all response instances
            fprintf('Computing %d response instances for a %d x %d cone mosaic with %d eye movements scanning a %d-frame stimulus sequence.\n', instancesNum, theConeMosaic.mosaicSize(1), theConeMosaic.mosaicSize(2), eyeMovementsNum, theOIsequence.length);
            tic
            [absorptionsCountSequence, absorptionsTimeAxis, photoCurrentSignals, photoCurrentTimeAxis] = ...
                    theConeMosaic.computeForOISequence(theOIsequence, ...
                    'emPaths', emPaths, ...
                    'currentFlag', true, ...
                    'newNoise', true ...
                    );
            fprintf('Response computation took %2.2f minutes\n', toc/60);
    
            %
            % Save responses for this stimStrengthIndex and chromaticDirectionIndex  here
            %
            
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
              absorptionsTimeAxis, ...
              absorptionsCountSequence(instancesToPlot,:,:,:), ...
              theOIsequence.oiTimeAxis, ...
              iL(idxL), iM(idxM), iS(idxS), coneStride, nominalStimulusStrengths(stimStrengthIndex));
         
          % Plot some L,M, and S cone response for the selected instances
          instancesToPlot = 1:1;
          coneStride = 20;  % how many cones to skip over
          plotAbsorptionsCountTimeSeries(...
              absorptionsTimeAxis, ...
              absorptionsCountSequence(instancesToPlot,:,:,:), ...
              theOIsequence.oiTimeAxis, ...
              iL, iM, iS,  coneStride, nominalStimulusStrengths(stimStrengthIndex));
          clear 'tmp'
        end % for contrastIndex
    end %  for chromaticDirectionIndex
    
end

function plotAbsorptionsCountTimeSeries(timeAxis, absorptions, stimTimeAxis, iL, iM, iS,  coneStride, stimulusStrength)
    
    instancesPlotted = size(absorptions,1);
    timePointsPlotted = size(absorptions,4);
    
    minAbsorptions = min(absorptions(:));
    maxAbsorptions = max(absorptions(:));
    
    % Extract cone responses by type
    iL = iL(1:coneStride:end);
    iM = iM(1:coneStride:end);
    iS = iS(1:coneStride:end);
    absorptions = permute(reshape(absorptions, [instancesPlotted size(absorptions,2)*size(absorptions,3) timePointsPlotted]), [1 3 2]);
    absorptionsByConeType {1} = reshape(permute(absorptions(:,:,iL), [1 3 2]), [instancesPlotted*numel(iL) timePointsPlotted]);
    absorptionsByConeType {2} = reshape(permute(absorptions(:,:,iM), [1 3 2]), [instancesPlotted*numel(iM) timePointsPlotted]);
    absorptionsByConeType {3} = reshape(permute(absorptions(:,:,iS), [1 3 2]), [instancesPlotted*numel(iS) size(absorptions,2)]);

    colors = [1 0 0; 0 1.0 0; 0 0.8 1];
    plotBackgroundColor = [0.1 0.1 0.1];
    barOpacity = 0.1;
    
    hFig = figure(); clf;
    set(hFig, 'Color', [0 0 0], 'Position', [10 10 650 900]);
    subplot('Position', [0.08 0.08 0.90 0.83]);
    
    hold on
    % Identify stimulus presentation times
    for k = 1:numel(stimTimeAxis)
        plot(stimTimeAxis(k)*[1 1]*1000, [minAbsorptions maxAbsorptions], 'k-', 'Color', [0.5 0.5 0.5]);
    end
        
    % Plot absorptions
    dt = timeAxis(2)-timeAxis(1);
    for coneType = 1:3
        for tIndex = 1:numel(timeAxis)
            quantaAtThisTimeBin = squeeze(absorptionsByConeType{coneType}(:,tIndex));
            plot([timeAxis(tIndex) timeAxis(tIndex)+dt]*1000, [quantaAtThisTimeBin(:) quantaAtThisTimeBin(:)], '-', 'LineWidth', 1.5, 'Color', [colors(coneType,:) barOpacity]);
        end
    end
    box on;
    xlabel('time (ms)', 'FontSize', 14);
    ylabel(sprintf('# of absorptions per %2.1f ms',dt*1000), 'FontSize', 14);
    set(gca, 'XLim', [timeAxis(1) timeAxis(end)+dt]*1000, 'YLim', [minAbsorptions maxAbsorptions], ...
             'Color', plotBackgroundColor, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], ...
             'FontSize', 12);
    title(sprintf('<Stimulus strength: %2.3f> Instances plotted: %d, Cones plotted: %d L-, %d M-, %d S-cone', stimulusStrength, instancesPlotted, numel(iL), numel(iM), numel(iS)), 'Color', [1 1 1], 'FontSize',14);
    drawnow;
end


function [theConeMosaic, eyeMovementsNum] = coneMosaicGenerate(mosaicParams, theOIsequence)
    % Default human mosaic
    theConeMosaic = coneMosaic;
    
    % Adjust size
    if isnan(mosaicParams.size)
        % Generate a human cone mosaic with 1L, 1M and 1S cone
        theConeMosaic.rows = 1;
        theConeMosaic.cols = 3;
        theConeMosaic.pattern = [2 3 4];
    else
        theConeMosaic.setSizeToFOV(mosaicParams.size);
    end
    
    % Set the noise
    theConeMosaic.noiseFlag = mosaicParams.photonNoise;

    % Set the integrationTime
    theConeMosaic.integrationTime = mosaicParams.integrationTime;
    
    % Generate the outer-segment object to be used by the coneMosaic
    theOuterSegment = osLinear();
    theOuterSegment.noiseFlag = mosaicParams.osNoise;
    
    % Set a custom timeStep, for @osLinear we do not need the default 0.1 msec
    theOuterSegment.timeStep = mosaicParams.osTimeStep;

    % Couple the outersegment object to the cone mosaic object
    theConeMosaic.os = theOuterSegment;

    % Generate eye movement sequence for all oi's
    stimulusSamplingInterval = theOIsequence.oiTimeAxis(2)-theOIsequence.oiTimeAxis(1);
    eyeMovementsNumPerOpticalImage = stimulusSamplingInterval/theConeMosaic.integrationTime;
    eyeMovementsNum = round(eyeMovementsNumPerOpticalImage*theOIsequence.length);
    
    if (eyeMovementsNum < 1)
        error('Less than 1 eye movement!!! \nStimulus sampling interval:%g Cone mosaic integration time: %g\n', stimulusSamplingInterval, theConeMosaic.integrationTime);
    else 
        fprintf('Optical image sequence contains %2.0f eye movements (%2.2f eye movements/oi)\n', eyeMovementsNum, eyeMovementsNumPerOpticalImage);
    end
end


function theOI = oiGenerate(oiParams)
    % Generate optics
    if ((isfield(oiParams, 'noOptics')) && (oiParams.noOptics == true))
        theOI = oiCreate('diffraction limited');
        optics = oiGet(theOI,'optics');           
        optics = opticsSet(optics,'fnumber',0);
        optics = opticsSet(optics, 'off axis method', 'skip');
        theOI = oiSet(theOI,'optics', optics);
    else
        theOI = oiCreate('wvf human');
    end
    
    % Remove Lens
    if ((isfield(oiParams, 'noLens')) && (oiParams.noLens == true))
        lens = oiGet(theOI,'lens');
        lens.density = 0;
        theOI = oiSet(theOI,'lens',lens);
    end
    
    % Set custom FOV
    if (isfield(oiParams, 'fieldOfViewDegs'))
        theOI = oiSet(theOI,'h fov',oiParams.fieldOfViewDegs);
    end
    
    % Set custom pupil diameter
    if (isfield(oiParams, 'pupilDiamMM'))
        focalLength = oiGet(theOI,'distance');
        desiredFNumber = focalLength/(oiParams.pupilDiamMm/1000);
        theOI  = oiSet(theOI ,'optics fnumber',desiredFNumber);
        pupilDiamMmCheck = 1000*oiGet(theOI,'optics aperture diameter');
        if (max(abs(pupilDiamMmCheck - oiParams.pupilDiamMm)) > 1e-8)
            error('Failed to set pupil diameter as expected');
        end
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