function generateSpatioTemporalPhotocurrentResponse()

    recomputeXTresponses = ~true;
    if (recomputeXTresponses)
        load('/Users/nicolas/Documents/MATLAB/projects/IBIOColorDetect/LconeProfile_null.mat')
        adaptationPhotonRates = signal/(5/1000);
        [C, ia, ic] = unique(coneXpos,'sorted');

        load('/Users/nicolas/Documents/MATLAB/projects/IBIOColorDetect/LconeProfile_stim.mat')
        pulsePhotonRates = signal/(5/1000) - adaptationPhotonRates;

        coneXpos = coneXpos(ia);
        adaptationPhotonRates = adaptationPhotonRates(ia);
        pulsePhotonRates = pulsePhotonRates(ia);
        meanAdaptationPhotonRates = meanActivation/(5/1000)*ones(size(adaptationPhotonRates));

        [timeAxis, XTresponse, XTnoisyResponses, XTresponseModulation, XTresponseSNR] = computeXTresponse(coneXpos, adaptationPhotonRates, pulsePhotonRates);
        [timeAxis, XTresponseMeanAdaptation, XTnoisyResponsesMeanAdaptation, XTresponseModulationMeanAdaptation, XTresponseMeanAdaptationSNR] = computeXTresponse(coneXpos, meanAdaptationPhotonRates, pulsePhotonRates);
    
        save('data.mat', 'timeAxis', 'coneXpos', 'XTresponse', 'XTresponseModulation', 'XTresponseSNR', 'XTnoisyResponses', 'XTresponseMeanAdaptation', 'XTresponseModulationMeanAdaptation', 'XTresponseMeanAdaptationSNR', 'XTnoisyResponsesMeanAdaptation');
    else
        
        load('data.mat', 'timeAxis', 'coneXpos', 'XTresponse', 'XTresponseModulation', 'XTresponseSNR', 'XTnoisyResponses', 'XTresponseMeanAdaptation', 'XTresponseModulationMeanAdaptation', 'XTresponseMeanAdaptationSNR', 'XTnoisyResponsesMeanAdaptation');
    end
    
    figure(2); clf;
    targetXPos = 0:0.05:2;
    for k = 1:numel(targetXPos)
        [~,coneIndex] = min(abs(targetXPos(k) - abs(coneXpos)));
        plot(timeAxis*1000, XTresponse(coneIndex,:), 'k-'); hold on;
        set(gca, 'XLim', [0 0.5]*1000);
    end

    
    
    
    figure(3); clf;
    subplot(1,3,1);
    contourLevelsDeltaMicroAmps = cat(2, [-15:0.5:-0.5], [0.5:0.5:15]);

    coneXposHiRes = coneXpos(1):0.002:coneXpos(end);
    timeAxisHiRes = timeAxis(1):1/1000:timeAxis(end);
    
    
    [X,Y] = meshgrid(coneXpos,timeAxis);
    [X2,Y2] = meshgrid(coneXposHiRes,timeAxisHiRes);
    
    XTresponseModulation = interp2(X,Y,XTresponseModulation',X2,Y2, 'spline',0);
    XTresponseModulationMeanAdaptation = interp2(X,Y,XTresponseModulationMeanAdaptation',X2,Y2, 'spline',0);
    
    contourf(coneXposHiRes, timeAxisHiRes*1000, XTresponseModulation, contourLevelsDeltaMicroAmps, 'LineColor', 'none');
    set(gca, 'CLim', [-10 10]);
    set(gca, 'YLim', [0 0.5]*1000);
    title('ecc-based adaptation');
    colorbar
    axis 'xy'
    
    subplot(1,3,2);
    contourf(coneXposHiRes, timeAxisHiRes*1000, XTresponseModulationMeanAdaptation, contourLevelsDeltaMicroAmps, 'LineColor', 'none');
    set(gca, 'CLim', [-10 10]);
    set(gca, 'YLim', [0 0.5]*1000);
    title('mean mosaic adaptation');
    colorbar
    axis 'xy'
    
    subplot(1,3,3);
    contourf(coneXposHiRes, timeAxisHiRes*1000, (XTresponseModulation-XTresponseModulationMeanAdaptation), contourLevelsDeltaMicroAmps, 'LineColor', 'none');
    set(gca, 'CLim', [-3 3]);
    set(gca, 'YLim', [0 0.5]*1000);
    title('ecc-based - mean mosaic adaptation')
    colorbar
    axis 'xy'
    
    cMap = brewermap(1024, '*RdBu');
    colormap(cMap);
    
    contourLevelsSNR = 0: 0.05 :5;
    
    
    XTresponseSNR = interp2(X,Y,XTresponseSNR',X2,Y2, 'spline',0);
    XTresponseMeanAdaptationSNR = interp2(X,Y,XTresponseMeanAdaptationSNR',X2,Y2, 'spline',0);
    
    figure(4); clf;
    subplot(1,3,1);
    contourf(coneXposHiRes, timeAxisHiRes*1000, XTresponseSNR, contourLevelsSNR,  'LineColor', 'none');
    set(gca, 'CLim', [0.0 3]);
    set(gca, 'YLim', [0 0.3]*1000);
    title('ecc-based adaptation SNR');
    set(gca, 'FontSize', 14);
    ylabel('\it time (msec)'); xlabel('\it space (degs)');
    colorbar
    axis 'xy'
    
    subplot(1,3,2);
    contourf(coneXposHiRes, timeAxisHiRes*1000, XTresponseMeanAdaptationSNR, contourLevelsSNR, 'LineColor', 'none');
    set(gca, 'CLim', [0.0 3]);
    set(gca, 'YLim', [0 0.3]*1000);
    title('mean mosaic adaptation SNR');
    xlabel('\it space (degs)');
    set(gca, 'FontSize', 14);
    colorbar
    axis 'xy'
    
    subplot(1,3,3);
    contourf(coneXposHiRes, timeAxisHiRes*1000, (XTresponseSNR./XTresponseMeanAdaptationSNR), contourLevelsSNR,  'LineColor', 'none');
    set(gca, 'CLim', [1/1.5 1.5]);
    set(gca, 'YLim', [0 0.3]*1000);
    title('ecc-based / mean mosaic adaptation SNR')
    xlabel('\it space (degs)');
    set(gca, 'FontSize', 14);
    colorbar
    axis 'xy'
    
    cMap = brewermap(1024, '*spectral');
    colormap(cMap);
    
end
    
function [timeAxis, XTresponse, XTnoisyResponses, XTresponseModulation, XTresponseSNR] = computeXTresponse(coneXpos, adaptationPhotonRates, pulsePhotonRates)
    figure(1); clf;
    plot(coneXpos, adaptationPhotonRates, 'k-'); hold on;
    plot(coneXpos, pulsePhotonRates, 'rs-')
    
    for xIndex = 1:numel(adaptationPhotonRates)
        fprintf('Running model for cone %d of %d\n', xIndex, numel(adaptationPhotonRates));
        adaptationPhotonRate = adaptationPhotonRates(xIndex);  % R*/cone/sec
        pulsePhotonRate = pulsePhotonRates(xIndex); % R*/cone/sec
        pulseDurationSeconds = 100/1000;
        photonsDeliveredDuringPulse = pulsePhotonRate * pulseDurationSeconds;

        cParams.spontaneousIsomerizationRate = 0; %  R*/c/s
        cParams.eccentricity  = 'foveal';
        cParams.useDefaultImplementation = true;
        cParams.noisyInstancesNum = 100;

        stimParams = struct(...
                'type', 'pulse', ...                                            % type of stimulus
                'timeSampleSeconds', 0.1/1000, ...                              % simulation time step
                'pulseDurationSeconds', pulseDurationSeconds, ...       % pulse duration in seconds
                'totalDurationSeconds', 1.0, ... % total duration of the stimulus
                'adaptationPhotonRate',  adaptationPhotonRate, ...      %  R*/c/s
                'photonsDeliveredDuringPulse', photonsDeliveredDuringPulse  ... %  R*/c during  the total pulse duration
            );

        % Design stimulus
        stimulus = designStimulus(stimParams, cParams.spontaneousIsomerizationRate);

         % The photonIntegrationTime affects the SNR of the cone
         % excitations. The longer it is the higher the SNR. 
         % Here we setting it to the pulse duration.
         vParams.photonIntegrationTime = pulseDurationSeconds;


        % Run model
        modelResponse = photocurrentModel(stimulus, cParams.eccentricity, cParams.noisyInstancesNum, ...
            vParams.photonIntegrationTime, cParams.useDefaultImplementation);
        
        % Differential response
        if (xIndex == 1)
            XTresponse = zeros(numel(adaptationPhotonRates), numel(modelResponse.membraneCurrent));
            XTnoisyResponses = zeros(numel(adaptationPhotonRates), size(modelResponse.noisyMembraneCurrents,1), size(modelResponse.noisyMembraneCurrents, 2));
            membraneNoiseSigma = zeros(1, numel(adaptationPhotonRates));
        end
        
        XTresponse(xIndex, :) = modelResponse.membraneCurrent;
        XTnoisyResponses(xIndex, :,:) = modelResponse.noisyMembraneCurrents;
        XTresponseModulation(xIndex, :) = modelResponse.membraneCurrent - modelResponse.membraneCurrent(1);
        membraneNoiseSigma(xIndex) = std(modelResponse.noisyMembraneCurrents(:), 0, 1);
    end
    
    XTresponseSNR = abs(XTresponseModulation)/mean(membraneNoiseSigma);
    
    timeAxis = modelResponse.timeAxis;
    
end

