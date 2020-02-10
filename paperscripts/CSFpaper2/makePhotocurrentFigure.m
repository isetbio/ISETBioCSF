function makePhotocurrentFigure
 
    % Define background luminace levels examined & the pedestal step (cd/m2)
    backgroundLuminances = [34];
    lumStep = 20;
    
    stimulusSamplingInterval = 20 / 1000;
    oiTimeAxis = -0.5:stimulusSamplingInterval:1;
    stimDuration = 500/1000;
    
    mosaicSize = nan;
    integrationTime = 0.1 / 1000;
    photonNoise = 'none';
    osNoise = 'none';
    osTimeStep = 0.1 / 1000;
    theConeMosaic = coneMosaicGenerate(...
        mosaicSize, photonNoise, osNoise, integrationTime, osTimeStep);
    
    theOI = oiCreate('human');
    
    modulationFunction = cell(numel(backgroundLuminances), 1);
    photocurrents = cell(numel(backgroundLuminances), 1);
    isomerizations = cell(numel(backgroundLuminances), 1);
    osLinearFilters = cell(numel(backgroundLuminances), 1);

    for iLum = 1:numel(backgroundLuminances)
        % Compute scene
        FOV = 1;
        theScene = uniformFieldSceneCreate(FOV, backgroundLuminances(iLum));
    
        oiBackground = oiCompute(theOI, theScene);
        oiModulated = oiBackground;

        modulationFunction{iLum} = zeros(1, numel(oiTimeAxis));
        modulationFunction{iLum}(oiTimeAxis >= 0.0 & oiTimeAxis < stimDuration) = ...
            lumStep / backgroundLuminances(iLum);
        theOIsequence = oiSequence(oiBackground, oiModulated, oiTimeAxis, ...
            modulationFunction{iLum});
    
        eyeMovementsNum = ...
            theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
        instancesNum = 1;
        emPaths = zeros(instancesNum, eyeMovementsNum, 2);
    
    end
    
    theOIsequence.visualize('montage');
    
    [isomerizations{iLum}, photocurrents{iLum}, ...
        osLinearFilters{iLum}] = ...
        theConeMosaic.computeForOISequence(theOIsequence, ...
        'emPaths', emPaths, 'currentFlag', true);
    
    timeAxis = theConeMosaic.timeAxis + theOIsequence.timeAxis(1);
    instanceNum = 1;
    figure(2); clf;
    subplot(1,2,1);
    plot(timeAxis, squeeze(isomerizations{iLum}(instanceNum,1,1,:)), 'r-');
    hold on;
    plot(timeAxis, squeeze(isomerizations{iLum}(instanceNum,1,2,:)), 'g-');
    plot(timeAxis, squeeze(isomerizations{iLum}(instanceNum,1,3,:)), 'b-');
    subplot(1,2,2);
    plot(timeAxis, squeeze(photocurrents{iLum}(instanceNum,1,1,:)), 'r-');
    hold on;
    plot(timeAxis, squeeze(photocurrents{iLum}(instanceNum,1,2,:)), 'g-');
    plot(timeAxis, squeeze(photocurrents{iLum}(instanceNum,1,3,:)), 'b-');
    set(gca, 'XLim', [0 1], 'YLim', [-86 -70]);

end

function theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, ...
    osNoise, integrationTime, osTimeStep)

    % Default human mosaic
    theConeMosaic = coneMosaic;

    % Adjust size
    if isnan(mosaicSize)
        % Generate a human cone mosaic with 1L, 1M and 1S cone
        theConeMosaic.rows = 1;
        theConeMosaic.cols = 3;
        theConeMosaic.pattern = [2 3 4];
    else
        theConeMosaic.setSizeToFOV(mosaicSize);
    end

    % Set the noise
    theConeMosaic.noiseFlag = photonNoise;

    % Set the integrationTime
    theConeMosaic.integrationTime = integrationTime;

    % Generate the outer-segment object to be used by the coneMosaic
    theOuterSegment = osLinear(); % osBioPhys();
    theOuterSegment.noiseFlag = osNoise;

    % Set a custom timeStep, for @osLinear we do not need the default 0.1 msec
    theOuterSegment.timeStep = osTimeStep;

    % Couple the outersegment object to the cone mosaic object
    theConeMosaic.os = theOuterSegment;

end

function uniformScene = uniformFieldSceneCreate(FOV, meanLuminance)
    uniformScene = sceneCreate('uniform equal photon', 128);

    % square scene with desired FOV
    uniformScene = sceneSet(uniformScene, 'wAngular', FOV);

    % 1 meter away
    uniformScene = sceneSet(uniformScene, 'distance', 1.0);

    % adjust radiance according to desired  mean luminance
    uniformScene = sceneAdjustLuminance(uniformScene, meanLuminance);
end
