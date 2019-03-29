function [theSNR, noiseEstimationLatency, peakEstimationLatency, modulationPeak, noiseSigma, responseAtPeakSigma, responsePeak] = computeSNR(timeAxis, noisyResponses, meanResponse)
    % measure noise properties from the last 0.5 seconds of the signal
    noiseEstimationLatency = timeAxis(end)-0.5;
    
    % Estimate sigma of noise
    tBinsForNoiseEstimation = find(timeAxis>=noiseEstimationLatency);
    noisyResponses2 = noisyResponses(:,tBinsForNoiseEstimation);
    noiseSigma = std(noisyResponses2(:),  0, 1);
    
    % Estimate mean of noise 
    noiseMean = mean(noisyResponses2(:));

    % Estimate mean of signal
    tBinsForSignalEstimation = find(timeAxis < noiseEstimationLatency);
    meanResponse = meanResponse(tBinsForSignalEstimation);
    
    % Compute modulation peak
    [modulationPeak, tPeakBinIndex] = max(abs(meanResponse-noiseMean));
    responsePeak = meanResponse(tPeakBinIndex);
    
    % Estimate sigma of response at peak
    peakEstimationLatency = timeAxis(tPeakBinIndex);
    noisyResponses3 = noisyResponses(:,tPeakBinIndex);
    responseAtPeakSigma = std(noisyResponses3(:),  0, 1);
    
    % Compute the SNR
    theSNR = modulationPeak / sqrt(0.5*(noiseSigma^2 + responseAtPeakSigma^2));
    
    % Show SNR components
    showSNRcomponents = false;
    if (showSNRcomponents)
        figure(1234); clf;
        plot(timeAxis(tPeakBinIndex)*[1 1], noiseMean+[0 modulationPeak], 'rs-', 'LineWidth', 1.5); hold on
        plot(timeAxis(tPeakBinIndex)*[1 1], meanResponse(tPeakBinIndex) + responseAtPeakSigma*[-1 1], 'gs-', 'LineWidth', 1.5);
        plot(1.2*[1 1], noiseMean+noiseSigma*[-1 1], 'co-', 'LineWidth', 1.5);
        plot(timeAxis, noisyResponses(1:20,:), 'k-'); 
        plot(timeAxis(tPeakBinIndex)*[1 1], noiseMean+[0 modulationPeak], 'rs-', 'LineWidth', 1.5);
        plot(timeAxis(tPeakBinIndex)*[1 1], meanResponse(tPeakBinIndex) + responseAtPeakSigma*[-1 1], 'gs-', 'LineWidth', 1.5);
        plot(1.0*[1 1], noiseMean+noiseSigma*[-1 1], 'co-', 'LineWidth', 1.5);
        legend({'modulationPeak', 'sigmaR', 'sigmaNoise'});
        drawnow;
    end
    
end
