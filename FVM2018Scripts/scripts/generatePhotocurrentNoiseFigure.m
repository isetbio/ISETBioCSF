function generatePhotocurrentNoiseFigure

FOV = 1;
meanLuminance = 30;
uniformScene = uniformFieldSceneCreate(FOV, meanLuminance);

theOI = oiCreate('human');
theOI = oiCompute(theOI, uniformScene);

integrationTime = 1/1000;
theMosaic = coneMosaicGenerate(nan, integrationTime);

fixationalEyeMovementsNum = 1.0/integrationTime;
nTrials = 500;
theEMPaths = generateTheEMPath(theMosaic, fixationalEyeMovementsNum, nTrials);
   


[theIsomerizations, thePhotocurrents, LMSfilters, noiseFreePhotocurrents] = ...
        theMosaic.compute(theOI, 'currentFlag',true, 'emPath', theEMPaths);

nTimeBins = size(theIsomerizations,4);
photocurrentNoiseSpectrum = zeros(nTrials, 3, nTimeBins);
    
hFig1 = figure(1); clf;
hFig2 = figure(2); clf;
set(hFig2, 'Position', [10 10 540 740], 'Color', [1 1 1]);
for coneTypeIndex = 1:3
    
    switch (coneTypeIndex)
        case 1
            color = 'r';
        case 2
            color = 'g';
        case 3
            color = 'b';
    end
    
    photocurrentNoise(:,1,coneTypeIndex,:) = thePhotocurrents(:,1,coneTypeIndex,:) - noiseFreePhotocurrents(coneTypeIndex);
    
    for iTrial = 1:nTrials
        photocurrentNoiseSpectrum(iTrial,coneTypeIndex,:) = reshape(...
            (abs((2/fixationalEyeMovementsNum) * fft(squeeze(photocurrentNoise(iTrial,1,coneTypeIndex,:))))).^2, ...
            [1 1 nTimeBins]);
    end
    meanNoiseSpectrum = squeeze(mean(photocurrentNoiseSpectrum, 1));
    
    nyquistFreq = 1/(2*integrationTime);
    deltaFreq = nyquistFreq/(0.5*nTimeBins);
    frequencyHz = deltaFreq * (0:(nTimeBins/2-1));
    
    figure(hFig1)
    subplot(3,2,(coneTypeIndex-1)*2+1);
    isomerizationResponse = squeeze(mean(theIsomerizations(:,1,coneTypeIndex,:),1));
    isomerizationRate(coneTypeIndex) = mean(isomerizationResponse)/integrationTime;
    plot(theMosaic.timeAxis, isomerizationResponse, '.-', 'Color', color);
    xlabel('time (sec)');
    ylabel('R*/cone/time bin');
    set(gca, 'YLim', [0 20], 'XLim', [0 1],  'FontSize', 14);
    
    subplot(3,2,(coneTypeIndex-1)*2+2);
    plot(theMosaic.timeAxis, squeeze(mean(thePhotocurrents(:,1,coneTypeIndex,:),1))-noiseFreePhotocurrents(coneTypeIndex), '.-', 'Color', color, 'LineWidth', 2.0);
    xlabel('time (sec)');
    ylabel('pA');
    set(gca, 'YLim', [-10 10], 'XLim', [0 1],  'FontSize', 16);
    
    figure(hFig2);
    subplot('Position', [0.18 0.59 0.78 0.37]);
    hold on
    nBins = size(LMSfilters,1);
    timeAxis = theMosaic.timeAxis;
    timeAxis = timeAxis(1:nBins);
    plot(timeAxis, squeeze(LMSfilters(:,coneTypeIndex)), 'k-', 'Color', color, 'LineWidth', 3.0);
    set(gca, 'YLim', [-0.1 0.65], 'XLim', [0 0.4], 'FontSize', 18, 'LineWidth', 0.7);
    xTicks = 0:0.05:0.4;
    yTicks = -.1:0.1:1;
    set(gca, 'XTick', xTicks , 'XTickLabel', sprintf('%.2f\n', xTicks ), ...
             'YTick', yTicks, 'YTickLabel', sprintf('%.2f\n', yTicks));
    grid on; box on
    xlabel('time (sec)', 'FontWeight', 'bold');
    ylabel(sprintf('single photon-on-background\nresponse (pA)'), 'FontWeight', 'bold');
    if (coneTypeIndex == 3)
        legend({...
            sprintf('L-cone (isomerization rate: %2.0f R*/sec)', isomerizationRate(1)), ...
            sprintf('M-cone (isomerization rate: %2.0f R*/sec)', isomerizationRate(2)), ...
            sprintf('S-cone (isomerization rate: %2.0f R*/sec)', isomerizationRate(3)), ...
            });
    end
    title(sprintf('pCurrent impulse response (bkgnd: %2.1f cd/m^{2})', meanLuminance));
    
    subplot('Position', [0.18 0.07 0.78 0.37]);
    loglog(frequencyHz, squeeze(meanNoiseSpectrum(coneTypeIndex,1:numel(frequencyHz))), '-', 'Color', color, 'LineWidth', 3.0);
    hold on
    set(gca, 'XTick', [1 3 10 30 100 300], 'YTick', [1e-3 1e-2 1e-1 1e0 1e1], 'XLim', [1 1000], 'YLim', [1e-3 1e1*0.5]);
    grid on; box on;
    if (coneTypeIndex == 3)
        legend({'L-cone', 'M-cone', 'S-cone'})
    end
    xlabel('frequency (Hz)', 'FontWeight', 'bold');
    ylabel('power (pA^{2} Hz^{-1})', 'FontWeight', 'bold');
    set(gca, 'XLim', [1 500],  'FontSize', 18, 'LineWidth', 0.7);
    title(sprintf('pCurrent noise (bkgnd: %2.1f cd/m^{2})', meanLuminance));
end

NicePlot.exportFigToPDF(sprintf('Photocurrent_%2.1fCDM2.pdf', meanLuminance), hFig2, 300);

end


    
    
function theConeMosaic = coneMosaicGenerate(mosaicSize, integrationTime)
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

    % Set the integrationTime
    theConeMosaic.integrationTime = integrationTime;
end

function theEMPath = generateTheEMPath(theMosaic, eyeMovementsNum, nTrials)
    fixEMobj = fixationalEM();
    fixEMobj.computeForConeMosaic(theMosaic, eyeMovementsNum, ...
        'nTrials', nTrials, ...
        'rSeed', 857);
    theEMPath = fixEMobj.emPos;
    
    if (1==2)
    fixEMobj.generateEMandMosaicComboVideo(...
            fixEMobj, theMosaic, ...
            'visualizedFOVdegs', 0.8*theMosaic.fov(1), ...
            'showMovingMosaicOnSeparateSubFig', true, ...
            'displaycrosshairs', true, ...
            'videoFileName', fullfile('../updatedComponentFigs', 'fixationalEMVideo.mp4'));
    end
        
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
