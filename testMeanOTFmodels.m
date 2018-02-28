function testMeanOTFmodels

    wavelengthsListToCompute = 400:30:700;

    recomputeOTFdata = ~true;
    if (recomputeOTFdata)
        centeringWavelength = 550;
        targetPupilDiamMM = 3.0;
        wavefrontSpatialSamples = 261*2+1;
        psfRange = 6*[-1 1];
        otfRange = [0 100];
        showTranslation = false;
        d = computeAllOTFs(wavelengthsListToCompute,centeringWavelength,targetPupilDiamMM, wavefrontSpatialSamples, psfRange, otfRange, showTranslation);
    else
        fprintf('Loading data\n');
        d = loadOTFs();
    end
    
    
    
    % Extract 1D OTF slices
    [xSfCyclesDeg, meanZcoeffMTFSlicesX, meanSubjectMTFSlicesX, subjectMTFSlicesX, ...
        meanZcoeffMTFSlicesY, meanSubjectMTFSlicesY, subjectMTFSlicesY] = ...
        extractMTFslices(d.xSfCyclesDeg, d.meanZcoeffOTF, d.meanSubjectOTF, d.subjectOTF);
    
    targetWavelength = 550;
    
    fprintf('Computing PCA on PSFs\n');
    % compute PSF matching scores
    pcaProjectionSpaceDim = size(d.subjectPSF,1);
    [PCAprojections, meanZcoeffProjection] = computePCAprojections(d.subjectPSF, d.meanZcoeffPSF, pcaProjectionSpaceDim);

    plotWeights = true;
    psfWaveWeightingSigma = 50*100;
    psfScores = computePSFmatchScore(PCAprojections, meanZcoeffProjection, wavelengthsListToCompute, targetWavelength, psfWaveWeightingSigma, plotWeights);
    
    mtfBiases = [0 0.05 0.1 0.2];
    mtfWaveWeightingSigmas = [1 10 30 60 100 1000];
    
    [totalScoresAcrossMethods, mtfScoresAcrossMethods, methodNames] = computeTotalScoresAcrossMTFMethods(psfScores, mtfBiases, mtfWaveWeightingSigmas, subjectMTF, targetMTF, wavelengthsListToCompute, targetWavelength);
    
    plotTotalScoresAcrossMTFMethods(totalScoresAcrossMethods, mtfScoresAcrossMethods, psfScores, methodNames);
    
    pause
    
    %[98 132 55 99 172 198 186 66];
    
    % compute OTF matching scores
    mtfWaveWeightingSigma = input('Enter spectral weighting sigma (OTF):');
    mtfBias = input('Enter spectral weighting bias (OTF):');
    mtfScores = computeMTFmatchScore(subjectMTF, targetMTF, wavelengthsListToCompute, targetWavelength, mtfWaveWeightingSigma, mtfBias, plotWeights);
     
    plotScores(mtfScores, psfScores, mtfWaveWeightingSigma, mtfBias);
     
    while (1)   
        whichSubject = input('Subject id to display:');
        plotRankedSubjectData(whichSubject, wavelengthsListToCompute, xSfCyclesDeg, d.xMinutes, ...
                meanZcoeffMTFSlicesX, meanSubjectMTFSlicesX, subjectMTFSlicesX, ...
                meanZcoeffMTFSlicesY, meanSubjectMTFSlicesY, subjectMTFSlicesY, ...
                squeeze(d.subjectPSF(whichSubject,:,:,:)), d.meanZcoeffPSF);
    end
    
    %plot2DPSFAndOTF(wavelengthsListToCompute, d.xMinutes, d.yMinutes,  d.xSfCyclesDeg, d.ySfCyclesDeg, d.meanZcoeffPSF, d.meanZcoeffOTF, d.meanSubjectOTF);
end

function plotTotalScoresAcrossMTFMethods(totalScoresAcrossMethods, mtfScoresAcrossMethods,psfScores, methodNames)
    mediansMTF = median(mtfScoresAcrossMethods,1);
    medians = median(totalScoresAcrossMethods,1);
    mads = mad(totalScoresAcrossMethods,0,1);
    
    [~,idx] = sort(psfScores, 'descend');
    %[~,idx] = sort(medians, 'descend');
    
    
    methodsNum = size(totalScoresAcrossMethods,1);
    nSubjects = size(totalScoresAcrossMethods,2);
    methodsColors = brewermap(methodsNum, 'Spectral');
    subjectIndices = 1:nSubjects;
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 1, ...
           'heightMargin',   0.06, ...
           'widthMargin',    0.001, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.08, ...
           'topMargin',      0.001);
       
    hFig = figure(11); clf;
    set(hFig, 'Position', [1 1 2560 560], 'Color', [1 1 1]);
    subplot('Position', subplotPosVectors(1,1).v);
    hold on;
    for methodIndex = 1:methodsNum
        plot(subjectIndices, squeeze(totalScoresAcrossMethods(methodIndex,:)), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', squeeze(methodsColors(methodIndex,:)));
    end
    set(gca, 'XLim', [0 nSubjects+1], 'YLim', [0 1], 'YTick', 0:0.1:1.0, 'XTick', 1:1:200, 'FontSize', 10);
    set(gca,'TickLength',[0.002, 0.002])
    hL = legend(methodNames, 'Orientation', 'horizontal', 'Location', 'northoutside');
    set(hL, 'FontSize', 12);
    xtickangle(60)
    ylabel('total score', 'FontSize', 14);
    grid on; box on;
    

    subplot('Position', subplotPosVectors(2,1).v);
    hold on
    plot(subjectIndices, medians(idx),'ko', 'LineWidth', 1.5, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 10);
    plot(subjectIndices, psfScores(idx), 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.5 0.8 1], 'MarkerSize', 6);
    plot(subjectIndices, mediansMTF(idx), 'ro-', 'LineWidth', 1.5, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', 6);

    for k = 1:nSubjects
        plot([k k], mads(idx(k))*[-1 1]+medians(idx(k)),'k-', 'LineWidth', 1.5);
    end
    hL = legend({'total score', 'psf score', 'otf score'}, 'Location', 'northeast');
    set(hL, 'FontSize', 14);
    ticks = subjectIndices;
    set(gca, 'XLim', [0 nSubjects+1], 'YLim', [0 1], 'YTick', 0:0.1:1.0, 'XTick', ticks, 'XTickLabel', sprintf('%d\n',subjectIndices(idx)), 'FontSize', 10);
    set(gca,'TickLength',[0.002, 0.002])
    xtickangle(60)
    xlabel('subject id', 'FontSize', 14);
    ylabel('scores', 'FontSize', 14);
    grid on; box on;
    NicePlot.exportFigToPDF(sprintf('scoresAllMethods.pdf'), hFig, 300);
end

function [totalScoresAcrossMethods, mtfScoresAcrossMethods, methodNames] = computeTotalScoresAcrossMTFMethods(psfScores, mtfBiases, mtfWaveWeightingSigmas, subjectMTF, targetMTF, wavelengthsListToCompute, targetWavelength)

    nSubjects = size(subjectMTF,1);
    methodIndex = 0;
    totalScoresAcrossMethods = zeros(numel(mtfBiases)*numel(mtfWaveWeightingSigmas), nSubjects);
    mtfScoresAcrossMethods = totalScoresAcrossMethods;
    
    for mtfBiasIndex = 1:numel(mtfBiases)
        for mtfSigmaIndex = 1:numel(mtfWaveWeightingSigmas)
            mtfScores = computeMTFmatchScore(subjectMTF, targetMTF, wavelengthsListToCompute, targetWavelength, ...
                mtfWaveWeightingSigmas(mtfSigmaIndex), mtfBiases(mtfBiasIndex), false);
            
            totalScores = sqrt(psfScores.^2+mtfScores.^2)/sqrt(2.0);
            methodIndex = methodIndex+1;
            totalScoresAcrossMethods(methodIndex,:) = totalScores;
            mtfScoresAcrossMethods(methodIndex,:) = mtfScores; 
            methodNames{methodIndex} = sprintf('b:%0.2f,s:%2.0f', mtfBiases(mtfBiasIndex),  mtfWaveWeightingSigmas(mtfSigmaIndex));
        end
    end
end

function plotScores(otfScores, psfScores, otfWaveWeightingSigma, otfBias)
    hFig = figure(32); clf;
    set(hFig, 'Position', [10 10 820 820], 'Color', [1 1 1]);
    
    for k = 1:numel(otfScores)
        text(otfScores(k), psfScores(k), sprintf('%d', k), 'Color', 'k', 'FontSize', 12);
    end
    xlabel('otf score');
    ylabel('psf score');
    set(gca, 'XTick', 0:0.1:1, 'YTick', 0:0.1:1, 'XLim', [-0.05 1.05], 'YLim', [-0.05 1.05],'FontSize', 16);
    axis 'square'; axis 'xy'; grid 'on'; box 'on'
    NicePlot.exportFigToPDF(sprintf('scores_otfSigma_%2.2f_otfBias_%2.2f.pdf', otfWaveWeightingSigma, otfBias), hFig, 300);
end

function [subjectPCAprojections, meanZcoeffPCAProjection] = computePCAprojections(subjectPSFs, meanZcoeffPSF, projectionSpaceDim)

    subjectsNum = size(subjectPSFs,1);
    N = size(subjectPSFs,2)*size(subjectPSFs,3);
    nWaves = size(subjectPSFs,4);
    
    subjectPCAprojections = zeros(nWaves,subjectsNum, projectionSpaceDim);
    meanZcoeffPCAProjection = zeros(nWaves, projectionSpaceDim);
    
    for iW = 1:nWaves 
        waveSubjectPSFs = squeeze(subjectPSFs(:,:,:,iW));
        waveSubjectPSFs = reshape(waveSubjectPSFs, [subjectsNum N]);
        
        waveMeanZcoeffPSF = squeeze(meanZcoeffPSF(:,:,iW));
        waveMeanZcoeffPSF = reshape(waveMeanZcoeffPSF, [1 N]);
        
        [~, ~, v] = svd(waveSubjectPSFs, 'econ');
%         reconstructedPSFs = u * s *v';
%         varianceExplained = (diag(s)).^2;
%         varianceExplained = varianceExplained / sum(varianceExplained);
%         explained = cumsum(varianceExplained);
%         explained = explained / max(explained)*100;

        
        for subjectIndex = 1:subjectsNum
           subjectPCAprojections(iW,subjectIndex,:) = squeeze(waveSubjectPSFs(subjectIndex,:) * v(:,1:projectionSpaceDim));
        end
        
        meanZcoeffPCAProjection(iW,:) = waveMeanZcoeffPSF * v(:,1:projectionSpaceDim);
    end
end

function psfScores = computePSFmatchScore(PCAprojections, meanZcoeffProjection, wavelengths, targetWavelength, otfWaveWeightingSigma, plotWeights)

    bias = 1;
    spectralWeighting = bias + (1-exp(-0.5*((wavelengths-targetWavelength)/otfWaveWeightingSigma).^2));
    spectralWeighting = spectralWeighting / max(spectralWeighting);
    
    if (plotWeights)
        hFig = figure(3); clf
        set(hFig, 'Position', [10 10 700 340], 'Color', [1 1 1]);
        subplot(1,2,1);
        stem(wavelengths, spectralWeighting, 'filled');
        set(gca, 'YLim', [-0.05 1.05], 'FontSize', 14);
        set(gca, 'XLim', [390 710]);
        axis 'square'; box on; grid on;
        xlabel('wavelength (nm)');
        ylabel('PSF weighting');
        drawnow
    end
    
    nWaves = size(meanZcoeffProjection,1);
    for iW = 1:nWaves 
        targetCoeffs = squeeze(meanZcoeffProjection(iW,:));
        subjectCoeffs = squeeze(PCAprojections(iW,:,:));
        distances(iW,:) = spectralWeighting(iW) * sqrt(sum((bsxfun(@minus, subjectCoeffs, targetCoeffs)).^2,2));
    end
    weightedMeanDistances = squeeze(sum(distances,1));
    
    psfRMS = (weightedMeanDistances-min(weightedMeanDistances))/(max(weightedMeanDistances)-min(weightedMeanDistances));
    psfScores = 1-psfRMS;
end

function otfScores = computeMTFmatchScore(subjectMTFs, targetMTF, wavelengths, targetWavelength, otfWaveWeightingSigma, otfBias, plotWeights)

    spectralWeighting = otfBias + (1-otfBias)*exp(-0.5*((wavelengths-targetWavelength)/otfWaveWeightingSigma).^2);
    if (plotWeights)
        hFig = figure(3); 
        subplot(1,2,2)
        stem(wavelengths, spectralWeighting, 'filled');
        set(gca, 'YLim', [-0.05 1.05], 'FontSize', 14);
        set(gca, 'XLim', [390 710]);
        axis 'square'; box on; grid on;
        xlabel('wavelength (nm)');
        ylabel('OTF weighting');
        drawnow;
        NicePlot.exportFigToPDF(sprintf('weights_otfSigma_%2.2f_otfBias_%2.2f.pdf',otfWaveWeightingSigma, otfBias), hFig, 300);
    end
    
    nSubjects = size(subjectMTFs,1);
    otfRMS = zeros(1, nSubjects);

    for iSubject = 1:nSubjects
        thisSubjectMTF = squeeze(subjectMTFs(iSubject,:,:,:));
        residual = (thisSubjectMTF-targetMTF).^2;
        weightedResidual = bsxfun(@times, residual, reshape(spectralWeighting, [1 1 numel(spectralWeighting)]));
        otfRMS(iSubject) = sqrt(sum(weightedResidual(:)));
    end

    otfRMS = (otfRMS-min(otfRMS))/(max(otfRMS)-min(otfRMS));
    otfScores = 1-otfRMS;
end

function [xSfCyclesDeg, meanZcoeffOTFSlicesX, meanSubjectOTFSlicesX, subjectOTFSlicesX, ...
          meanZcoeffOTFSlicesY, meanSubjectOTFSlicesY, subjectOTFSlicesY] = ...
    extractMTFslices( xSfCyclesDeg, meanZcoeffOTF, meanSubjectOTF, subjectOTF)
    midRow = floor(size(subjectOTF,2)/2)+1;
    
    meanZcoeffOTFSlicesX = abs(squeeze(meanZcoeffOTF(midRow,midRow:end,:)));
    meanSubjectOTFSlicesX = abs(squeeze(meanSubjectOTF(midRow,midRow:end,:)));
    subjectOTFSlicesX = abs(squeeze(subjectOTF(:,midRow,midRow:end,:)));
    
    meanZcoeffOTFSlicesY = abs(squeeze(meanZcoeffOTF(midRow:end,midRow,:)));
    meanSubjectOTFSlicesY = abs(squeeze(meanSubjectOTF(midRow:end,midRow,:)));
    subjectOTFSlicesY = abs(squeeze(subjectOTF(:,midRow:end,midRow,:)));
    
    xSfCyclesDeg = xSfCyclesDeg(midRow:end);

end

function plotRankedSubjectData(subjectIndex, wavelengthsListToCompute, ...
    xSfCyclesDeg, xMinutes, meanZcoeffOTFSlicesX, meanSubjectOTFSlicesX, subjectOTFSlicesX, ...
    meanZcoeffOTFSlicesY, meanSubjectOTFSlicesY, subjectOTFSlicesY, ...
    rankedSubjectPSF, meanZcoeffPSF)

    nWaves = size(subjectOTFSlicesX,3);
    waveStride = 3;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 3, ...
           'colsNum', floor(nWaves/waveStride), ...
           'heightMargin',   0.03, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.01, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.01, ...
           'topMargin',      0.03);
       
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1670 600], 'Color', [1 1 1]);
    
    idx = find(wavelengthsListToCompute==550);
    targetMod = mod(idx-1,waveStride);
    
    iiw = 0;
    for iW = 1:nWaves 
        if (mod(iW-1,waveStride) ~= targetMod) || (iiw+1>size(subplotPosVectors,2))
            %fprintf('Skipping %dnm\n', wavelengthsListToCompute(iW));
            continue;
        end
        iiw = iiw + 1; 
        
        ySliceColor = [1 0 0];
        ySliceColor2 = [1 0.5 0.5];
        xSliceColor = [0. 0. 1.0];
        xSliceColor2 = [0.5 0.5 1.0];
        
        subplot('Position', subplotPosVectors(1,iiw).v);
        plot(xSfCyclesDeg, squeeze(meanSubjectOTFSlicesX(:,iW)), ':', 'Color', xSliceColor, 'LineWidth', 2.0);
        hold on;
        plot(xSfCyclesDeg, squeeze(meanSubjectOTFSlicesY(:,iW)), ':', 'Color', ySliceColor, 'LineWidth', 2.0);
        plot(xSfCyclesDeg, squeeze(subjectOTFSlicesX(subjectIndex,:,iW)), '-', 'Color', xSliceColor , 'LineWidth', 2.0);
        plot(xSfCyclesDeg, squeeze(subjectOTFSlicesY(subjectIndex,:,iW)), '-', 'Color', ySliceColor , 'LineWidth', 2.0);
        
        plot(xSfCyclesDeg, squeeze(subjectOTFSlicesX(1:2:end,:,iW)), '-', 'Color', [0.6 0.6 0.6]);
        plot(xSfCyclesDeg, squeeze(subjectOTFSlicesY(1:2:end,:,iW)), '-', 'Color', [0.6 0.6 0.6]);
        
        plot(xSfCyclesDeg, squeeze(subjectOTFSlicesX(subjectIndex,:,iW)), '-', 'Color', xSliceColor2 , 'LineWidth', 4.0);
        plot(xSfCyclesDeg, squeeze(subjectOTFSlicesY(subjectIndex,:,iW)), '-', 'Color', ySliceColor2 , 'LineWidth', 4.0);
        plot(xSfCyclesDeg, squeeze(meanSubjectOTFSlicesX(:,iW)), '-', 'Color', xSliceColor2, 'LineWidth', 4.0);
        plot(xSfCyclesDeg, squeeze(meanSubjectOTFSlicesY(:,iW)), '-', 'Color', ySliceColor2, 'LineWidth', 4.0);
        plot(xSfCyclesDeg, squeeze(meanSubjectOTFSlicesX(:,iW)), ':', 'Color', xSliceColor, 'LineWidth', 2.0);
        plot(xSfCyclesDeg, squeeze(meanSubjectOTFSlicesY(:,iW)), ':', 'Color', ySliceColor, 'LineWidth', 2.0);
        plot(xSfCyclesDeg, squeeze(subjectOTFSlicesX(subjectIndex,:,iW)), '-', 'Color', xSliceColor , 'LineWidth', 2.0);
        plot(xSfCyclesDeg, squeeze(subjectOTFSlicesY(subjectIndex,:,iW)), '-', 'Color', ySliceColor , 'LineWidth', 2.0);
        
        hold off
        set(gca, 'FontSize', 14, 'XTick', 0:20:100, 'XLim', [0 100], 'YLim', [0 1]);
        if (iW > 1)
            set(gca, 'YTickLabel', {});
            set(gca, 'XTickLabel', {});
        end
        axis 'square'
        grid on; box on;
        legend({'mean(x)', 'mean(y)', sprintf('#%d (x)', subjectIndex), sprintf('#%d (y)', subjectIndex)});
        title(sprintf('%d nm', wavelengthsListToCompute(iW)));
        
        
        subplot('Position', subplotPosVectors(2,iiw).v);
        thePSF = squeeze(rankedSubjectPSF(:,:,iW));
        imagesc(xMinutes, xMinutes, thePSF);
        hold on;
        plot([xMinutes(1) xMinutes(end)], [0 0], 'g-');
        plot([0 0], [xMinutes(1) xMinutes(end)],  'g-');
        hold off;
        axis 'image'; 
        if (iW > 1)
            set(gca, 'YTickLabel', {});
        end
        set(gca, 'XTickLabel', {});
        set(gca, 'FontSize', 12);
        title(sprintf('#%d PSF',subjectIndex));
        
        subplot('Position', subplotPosVectors(3,iiw).v);
        thePSF = squeeze(meanZcoeffPSF(:,:,iW));
        imagesc(xMinutes, xMinutes, thePSF);
        hold on;
        plot([xMinutes(1) xMinutes(end)], [0 0], 'g-');
        plot([0 0], [xMinutes(1) xMinutes(end)],  'g-');
        hold off;
        axis 'image'; 
        if (iW > 1)
            set(gca, 'YTickLabel', {});
        end
        set(gca, 'XTickLabel', {});
        set(gca, 'FontSize', 12);
        title('meanZ PSF');
        
    end
    colormap(bone(1024));
    
    NicePlot.exportFigToPNG(sprintf('subject%d.png', subjectIndex),hFig, 300);
end

function plot2DPSFAndOTF(wavelengthsListToCompute, xMinutes, yMinutes, xSfCyclesDeg, ySfCyclesDeg, thePSF, theOTF, theMeanSubjectOTF)
    figure(2); clf;
    nWaves = numel(wavelengthsListToCompute);
    
    for iW = 1:nWaves  
        subplot(3, nWaves,iW);
        imagesc(xMinutes, yMinutes, squeeze(thePSF(:,:,iW)));
        title(sprintf('%d nm', wavelengthsListToCompute(iW)));
        axis 'image'; xlabel('minutes');
        if (iW > 1)
            set(gca, 'YTickLabel', {});
        end
        set(gca, 'FontSize', 12);
        
        subplot(3, nWaves,iW+nWaves);
        imagesc(xSfCyclesDeg, ySfCyclesDeg, squeeze(theOTF(:,:,iW)));
        axis 'image'; xlabel('c/deg');
        if (iW > 1)
            set(gca, 'YTickLabel', {});
        end
        set(gca, 'FontSize', 12);
        
        subplot(3, nWaves,iW+2*nWaves);
        imagesc(xSfCyclesDeg, ySfCyclesDeg, squeeze(theMeanSubjectOTF(:,:,iW)));
        axis 'image'; xlabel('c/deg');
        if (iW > 1)
            set(gca, 'YTickLabel', {});
        end
        set(gca, 'FontSize', 12);
        
    end
    colormap(gray(1024));
end

function d = loadOTFs()
    load('allPSFdata.mat', ...
        'meanZcoeffPSF', 'meanZcoeffOTF',  ...
        'meanSubjectOTF',  'subjectPSF', 'subjectOTF', ...
        'xMinutes', 'yMinutes', 'xSfCyclesDeg', 'ySfCyclesDeg');
    

    d.xMinutes = xMinutes;
    d.yMinutes = yMinutes;
    d.xSfCyclesDeg = xSfCyclesDeg;
    d.ySfCyclesDeg = ySfCyclesDeg;
    d.meanZcoeffPSF = meanZcoeffPSF;
    d.meanZcoeffOTF = meanZcoeffOTF;
    d.meanSubjectOTF = meanSubjectOTF;
    d.subjectPSF = subjectPSF;
    d.subjectOTF = subjectOTF;

end

function d = computeAllOTFs(wavelengthsListToCompute,centeringWavelength,targetPupilDiamMM, wavefrontSpatialSamples, psfRange, otfRange, showTranslation)
    % Load the ZernikeCoeffs for the target pupil size
    [Zcoeffs_SampleMean, Zcoeffs_S, subject_coeffs] = wvfLoadThibosVirtualEyes(targetPupilDiamMM);
    ZcoeffSubjects = subject_coeffs.bothEyes;
    subjectsNum = size(ZcoeffSubjects,2);
    
    % Compute the meanZcoeff PSF and OTF
    fprintf('computing mean Zcoeff PSF/OTF\n');
    [meanZcoeffPSF, meanZcoeffOTF, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes] = ...
        computePSFandOTF(Zcoeffs_SampleMean, wavelengthsListToCompute, wavefrontSpatialSamples, targetPupilDiamMM, centeringWavelength, showTranslation);
    
    % Only keep PSF data within the psfRange
    psfColsToKeep = find(abs(xMinutes) < max(psfRange));
    psfRowsToKeep = find(abs(yMinutes) < max(psfRange));
    xMinutes = xMinutes(psfColsToKeep);
    yMinutes = yMinutes(psfRowsToKeep);
    meanZcoeffPSF = meanZcoeffPSF(psfRowsToKeep, psfColsToKeep,:);
    
    % Only keep OTF data within the otfRange
    otfColsToKeep = find(abs(xSfCyclesDeg) < max(otfRange));
    otfRowsToKeep = find(abs(ySfCyclesDeg) < max(otfRange));
    xSfCyclesDeg  = xSfCyclesDeg(otfColsToKeep);
    ySfCyclesDeg  = ySfCyclesDeg(otfRowsToKeep);
    meanZcoeffOTF = meanZcoeffOTF(otfRowsToKeep, otfColsToKeep,:);
    
    % Preallocate memory to store the subject data
    subjectPSF = zeros(subjectsNum, size(meanZcoeffPSF,1), size(meanZcoeffPSF, 2), size(meanZcoeffPSF, 3));
    subjectOTF = zeros(subjectsNum, size(meanZcoeffOTF,1), size(meanZcoeffOTF, 2), size(meanZcoeffOTF, 3));
    
    % Compute individual subject data
    parfor subjectIndex = 1:subjectsNum
        fprintf('computing PSF/OTF for subject %d\n', subjectIndex);
        z = squeeze(ZcoeffSubjects(:,subjectIndex));    
        [psf, otf, ~, ~, ~, ~] = ...
            computePSFandOTF(z, wavelengthsListToCompute, wavefrontSpatialSamples, targetPupilDiamMM, centeringWavelength, showTranslation);
        subjectPSF(subjectIndex,:,:,:) = psf(psfRowsToKeep, psfColsToKeep,:);
        subjectOTF(subjectIndex,:,:,:) = otf(otfRowsToKeep, otfColsToKeep,:);
    end
    % Compute mean across subjects OTF and MTF
    meanSubjectOTF = squeeze(mean(subjectOTF,1));
    
    
    % Save data
    save('allPSFdata.mat', ...
        'meanZcoeffPSF', 'meanZcoeffOTF',  ...
        'meanSubjectOTF', 'subjectPSF', 'subjectOTF', ...
        'xMinutes', 'yMinutes', 'xSfCyclesDeg', 'ySfCyclesDeg', ...
        '-v7.3');
    fprintf('data saved\n');
    
    d.xMinutes = xMinutes;
    d.yMinutes = yMinutes;
    d.xSfCyclesDeg = xSfCyclesDeg;
    d.ySfCyclesDeg = ySfCyclesDeg;
    d.meanZcoeffPSF = meanZcoeffPSF;
    d.meanZcoeffOTF = meanZcoeffOTF;
    d.meanSubjectOTF = meanSubjectOTF;
    d.subjectPSF = subjectPSF;
    d.subjectOTF = subjectOTF;
end