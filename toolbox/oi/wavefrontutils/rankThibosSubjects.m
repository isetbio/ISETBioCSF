function rankThibosSubjects()

    % Recompute or load previously computed OTFs/PSFs for all 200 Thibos subjects
    recomputeOTFdata = ~true;
    [d,wavelengthsListToCompute, focusWavelength] = generateMultiSpectralOTFs(recomputeOTFdata);
    
    % Generate spectral weights for the PSF residuals
    psfSpectralWeights = computePSFspectralWeighting(wavelengthsListToCompute, focusWavelength);
    
    % Compute the PSF matching scores
    psfScores = computePSFmatchScore(d.subjectProjections, d.meanZcoeffProjection, psfSpectralWeights);
    
    % Compute the spectral weights for the MTF residuals
    [mtfSpetralWeights, weightingMethodNames] = computeMTFspectralWeigting(wavelengthsListToCompute, focusWavelength);
    
    % Use the magnitudes for the OTF score (MTF)
    [mtfScoresMedianAcrossMethods, mtfScoresDifferentMethods] = computeMTFmatchScores(abs(d.subjectOTF), abs(d.meanSubjectOTF), mtfSpetralWeights);
    
    rankedScore = 'PSFscore';
    %rankedScore = 'medianMTFscore';
    %rankedScore = 'totalScore';
    
    switch (rankedScore)
        case 'PSFscore'
            % Plot the subjects scores ranked based on the psfScore
            plotRankedSubjects(psfScores, mtfScoresMedianAcrossMethods, mtfScoresDifferentMethods, weightingMethodNames, 'psf score', 'mtf score (mean)', 'b', 'r');
    	case 'medianMTFscore'
            % Plot the subjects scores ranked based on the median across MTF
            plotRankedSubjects(mtfScoresMedianAcrossMethods, psfScores, mtfScoresDifferentMethods, weightingMethodNames, 'mtf score (mean)', 'psf score', 'r', 'b');
        case 'totalScore'
            totalScore = sqrt(mtfScoresMedianAcrossMethods.^2 + psfScores.^2)/sqrt(2.0);
            % Plot the subjects scores ranked based on the median across MTF
            plotRankedSubjects(totalScore, psfScores, mtfScoresDifferentMethods, weightingMethodNames, 'total score', 'psf score', 'k', 'b');
    end
    
    % Plot individual subjects
    [xSfCyclesDeg, meanZcoeffMTFSlicesX, meanSubjectMTFSlicesX, subjectMTFSlicesX, ...
        meanZcoeffMTFSlicesY, meanSubjectMTFSlicesY, subjectMTFSlicesY] = ...
        extractMTFslices(d.xSfCyclesDeg, d.meanZcoeffOTF, d.meanSubjectOTF, d.subjectOTF);
    
    while (1)   
        whichSubject = input('Subject id to display:');
        plotRankedSubjectData(whichSubject, wavelengthsListToCompute, xSfCyclesDeg, d.xMinutes, ...
                meanZcoeffMTFSlicesX, meanSubjectMTFSlicesX, subjectMTFSlicesX, ...
                meanZcoeffMTFSlicesY, meanSubjectMTFSlicesY, subjectMTFSlicesY, ...
                squeeze(d.subjectPSF(whichSubject,:,:,:)), d.meanZcoeffPSF);
    end
    
    %plotSPectralWeights(wavelengthsListToCompute, psfSpectralWeights, mtfSpectralWeights);
end

% Method to extract 1D MTF slices
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

% Method to compute the MTF matching scores
function [mtfMedianScores, mtfScores] = computeMTFmatchScores(subjectMTFs, targetMTF, spectralWeightings)
    fprintf('\nComputing MTF scores ...');
    weightingsNum = size(spectralWeightings,1) * size(spectralWeightings,2);

    nSubjects = size(subjectMTFs,1);
    mtfRMS = zeros(weightingsNum, nSubjects);
    weightIndex = 0;
    
    % Compute RMS residuals for all spectral weighting methods
    for i = 1:size(spectralWeightings,1)
    for j = 1:size(spectralWeightings,2)
        weights = squeeze(spectralWeightings(i,j,:));
        weightIndex = weightIndex+1;
        for iSubject = 1:nSubjects
            thisSubjectMTF = squeeze(subjectMTFs(iSubject,:,:,:));
            residual = (thisSubjectMTF-targetMTF).^2;
            weightedResidual = bsxfun(@times, residual, reshape(weights, [1 1 numel(weights)]));
            mtfRMS(weightIndex,iSubject) = sqrt(sum(weightedResidual(:)));
        end
    end
    end
    
    % Median of the rms residuals across all methods (without separate normalization for different methods)
    mtfRMSmedian = median(mtfRMS,1);
    
    % Normalize the median of the RMS residuals to [0 1]
    mtfRMSmedian  = (mtfRMSmedian -min(mtfRMSmedian (:)))/(max(mtfRMSmedian (:))-min(mtfRMSmedian(:)));
    mtfMedianScores = 1-mtfRMSmedian;
    
    % Scores for each method separately
    mtfScores = zeros(weightingsNum, nSubjects);
    for weightIndex = 1:weightingsNum
        w = squeeze(mtfRMS(weightIndex,:));
        % Normalize rms separately
        w = (w-min(w(:)))/(max(w(:))-min(w(:)));
        mtfScores(weightIndex,:) = 1-w;
    end
    fprintf('Done !\n');
end

% Method to compute the PSF matching scores.
% PSF score is computed by weigted average of the 
function psfScores = computePSFmatchScore(PCAprojections, meanZcoeffProjection, spectralWeighting)
    fprintf('\nComputing PSF scores ...');
    nWaves = size(meanZcoeffProjection,1);
    for iW = 1:nWaves 
        targetCoeffs = squeeze(meanZcoeffProjection(iW,:));
        subjectCoeffs = squeeze(PCAprojections(iW,:,:));
        distances(iW,:) = spectralWeighting(iW) * sqrt(sum((bsxfun(@minus, subjectCoeffs, targetCoeffs)).^2,2));
    end
    weightedMeanDistances = squeeze(sum(distances,1));
    
    psfRMS = (weightedMeanDistances-min(weightedMeanDistances))/(max(weightedMeanDistances)-min(weightedMeanDistances));
    psfScores = 1-psfRMS;
    fprintf('Done !\n');
end

% Method to compute the spectral weights for the MTF residuals
function [spectralWeightings, weightingMethodNames] = computeMTFspectralWeigting(wavelengths, targetWavelength)
    biases = [0 0.05 0.1 0.2];
    sigmas = [1 10 30 60 100 1000];
    
    spectralWeightings = zeros(numel(biases), numel(sigmas), numel(wavelengths));
    weightingMethodNames = {};
    for biasIndex = 1:numel(biases)
        bias =  biases(biasIndex);
        for sigmaIndex = 1:numel(sigmas)
            sigma = sigmas(sigmaIndex);
            weightingMethodNames{numel(weightingMethodNames)+1} = sprintf('b:%2.2f, s:%2.0f', bias, sigma);
            w = bias + (1-bias)*exp(-0.5*((wavelengths-targetWavelength)/sigma).^2);
            spectralWeightings(biasIndex, sigmaIndex,:) = w/sum(w);
        end
    end
end

% Method to compute the spectral weights for the PSF residuals
function spectralWeighting = computePSFspectralWeighting(wavelengths, targetWavelength)
    bias = 1; sigma = 1000;
    spectralWeighting = bias + (1-exp(-0.5*((wavelengths-targetWavelength)/sigma).^2));
    spectralWeighting = spectralWeighting / sum(spectralWeighting);
end

% Method to compute the projections of all subjects PSFs onto the space spanned by their SVs
function [subjectSVDprojections, meanZcoeffSVDProjection] = computeSVDprojections(subjectPSFs, meanZcoeffPSF, projectionSpaceDim)

    subjectsNum = size(subjectPSFs,1);
    N = size(subjectPSFs,2)*size(subjectPSFs,3);
    nWaves = size(subjectPSFs,4);
    
    subjectSVDprojections = zeros(nWaves,subjectsNum, projectionSpaceDim);
    meanZcoeffSVDProjection = zeros(nWaves, projectionSpaceDim);
    
    for iW = 1:nWaves 
        waveSubjectPSFs = squeeze(subjectPSFs(:,:,:,iW));
        waveSubjectPSFs = reshape(waveSubjectPSFs, [subjectsNum N]);
        waveMeanZcoeffPSF = squeeze(meanZcoeffPSF(:,:,iW));
        waveMeanZcoeffPSF = reshape(waveMeanZcoeffPSF, [1 N]);
        
        % Compute SVD
        [~, ~, v] = svd(waveSubjectPSFs, 'econ');
%         reconstructedPSFs = u * s *v';
%         varianceExplained = (diag(s)).^2;
%         varianceExplained = varianceExplained / sum(varianceExplained);
%         explained = cumsum(varianceExplained);
%         explained = explained / max(explained)*100;

        for subjectIndex = 1:subjectsNum
           subjectSVDprojections(iW,subjectIndex,:) = squeeze(waveSubjectPSFs(subjectIndex,:) * v(:,1:projectionSpaceDim));
        end 
        meanZcoeffSVDProjection(iW,:) = waveMeanZcoeffPSF * v(:,1:projectionSpaceDim);
    end
end


% Method to compute (or load previously computed) multispectral OTF/PSFs
function [d, wavelengthsListToCompute, focusWavelength] = generateMultiSpectralOTFs(recomputeOTFdata)
    % List of wavelengths to compute
    wavelengthsListToCompute = 400:10:700;
    focusWavelength = 550;
    
    if (recomputeOTFdata)
        targetPupilDiamMM = 3.0;
        wavefrontSpatialSamples = 261*2+1;
        psfRange = 6*[-1 1];
        otfRange = [0 100];
        showTranslation = false;
        d = computeOTFsAndPSFsForAllThibosSubjects(wavelengthsListToCompute, focusWavelength, targetPupilDiamMM, wavefrontSpatialSamples, psfRange, otfRange, showTranslation);
    else
        d = loadPreComputedOTFs();
    end
end

% Method to compute and save the multi-spectra OTF and PSFs for all Thibos's subjects
function d = computeOTFsAndPSFsForAllThibosSubjects(wavelengthsListToCompute,centeringWavelength,targetPupilDiamMM, wavefrontSpatialSamples, psfRange, otfRange, showTranslation)
    % Load the ZernikeCoeffs for the target pupil size
    [Zcoeffs_SampleMean, ~, subject_coeffs] = wvfLoadThibosVirtualEyes(targetPupilDiamMM);
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
    
    fprintf('Computing PCA on PSFs\n');
    % Choose to project into the full 200-dimensional SVD space
    pcaProjectionSpaceDim = size(subjectPSF,1);
    [subjectProjections, meanZcoeffProjection] = computeSVDprojections(subjectPSF, meanZcoeffPSF, pcaProjectionSpaceDim);
        
    fprintf('\nSaving data ...');
    % Save data
    save('allPSFdata.mat', ...
        'meanZcoeffPSF', 'meanZcoeffOTF',  ...
        'meanSubjectOTF', 'subjectPSF', 'subjectOTF', ...
        'subjectProjections', 'meanZcoeffProjection', ...
        'xMinutes', 'yMinutes', 'xSfCyclesDeg', 'ySfCyclesDeg', ...
        '-v7.3');
    fprintf('Done!\n');
    
    d.xMinutes = xMinutes;
    d.yMinutes = yMinutes;
    d.xSfCyclesDeg = xSfCyclesDeg;
    d.ySfCyclesDeg = ySfCyclesDeg;
    d.meanZcoeffPSF = meanZcoeffPSF;
    d.meanZcoeffOTF = meanZcoeffOTF;
    d.meanSubjectOTF = meanSubjectOTF;
    d.subjectPSF = subjectPSF;
    d.subjectOTF = subjectOTF;
    d.subjectProjections = subjectProjections;
    d.meanZcoeffProjection = meanZcoeffProjection;
end

% Method to load the precomputed multi-spectra OTF and PSFs for all Thibos's subjects
function d = loadPreComputedOTFs()
    fprintf('\nLoading data ...');
    load('allPSFdata.mat', ...
        'meanZcoeffPSF', 'meanZcoeffOTF',  ...
        'meanSubjectOTF',  'subjectPSF', 'subjectOTF', ...
        'subjectProjections', 'meanZcoeffProjection', ...
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
    d.subjectProjections = subjectProjections;
    d.meanZcoeffProjection = meanZcoeffProjection;
    fprintf('Done !\n');
end


% Plotting methods
% Method to plot the ranked scores
function plotRankedSubjects(rankingScores, otherScores, mtfScoresDiffMethods, weightingMethodNames, rankingScoreLabel, otherScoreLabel, rankingScoreColor, otherScoreColor)
    fprintf('\nPlotting ranked subjects ...');
    [~,idx] = sort(rankingScores, 'descend');
    methodsNum = size(mtfScoresDiffMethods,1);
    nSubjects = size(mtfScoresDiffMethods,2);
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
    
    % The mtf scores for all examined wavelength weighting methods
    subplot('Position', subplotPosVectors(1,1).v);
    hold on;
    for methodIndex = 1:methodsNum
        plot(subjectIndices, squeeze(mtfScoresDiffMethods(methodIndex,:)), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', squeeze(methodsColors(methodIndex,:)));
    end
    set(gca, 'XLim', [0 nSubjects+1], 'YLim', [0 1], 'YTick', 0:0.1:1.0, 'XTick', 1:1:200, 'FontSize', 10);
    set(gca,'TickLength',[0.002, 0.002])
    hL = legend(weightingMethodNames, 'Orientation', 'horizontal', 'Location', 'northoutside');
    set(hL, 'FontSize', 12);
    xtickangle(60)
    ylabel('total score', 'FontSize', 14);
    grid on; box on;
     
    
    if (~strcmp(rankingScoreLabel, 'total score')) && (~strcmp(otherScoreLabel, 'total score'))
        % Compute the total score
        totalScores = sqrt(rankingScores.^2+otherScores.^2)/sqrt(2.0);
    else
        totalScores = [];
    end
    
    subplot('Position', subplotPosVectors(2,1).v);
    hold on
    switch (rankingScoreColor)
        case 'b'
            plot(subjectIndices, rankingScores(idx),'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.5 0.8 1], 'MarkerSize', 10);
        case 'r'
            plot(subjectIndices, rankingScores(idx),'ro-', 'LineWidth', 1.5, 'MarkerFaceColor', [1 0.5, 0.5], 'MarkerSize', 10);
        case 'k'
            plot(subjectIndices, rankingScores(idx),'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.5 0.5, 0.5], 'MarkerSize', 10);
        otherwise
            error('Unknown color for rankingScore')
    end
    
    switch (otherScoreColor)
        case 'r'
            plot(subjectIndices, otherScores(idx),  'ro-', 'LineWidth', 1.5, 'MarkerFaceColor', [1 0.5, 0.5], 'MarkerSize', 10);
        case 'b'
            plot(subjectIndices, otherScores(idx),  'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.5 0.8 1], 'MarkerSize', 10);
        case 'k'
            plot(subjectIndices, otherScores(idx),  'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.5 0.5, 0.5], 'MarkerSize', 10);
    end
            
    if (~isempty(totalScores))
        plot(subjectIndices, totalScores(idx),  'ko', 'LineWidth', 1.5, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8);
        hL = legend({rankingScoreLabel, otherScoreLabel, 'total score'}, 'Location', 'northeast');
    else
        hL = legend({rankingScoreLabel, otherScoreLabel}, 'Location', 'northeast');
    end
    
    
    set(hL, 'FontSize', 14);
    ticks = subjectIndices;
    set(gca, 'XLim', [0 nSubjects+1], 'YLim', [0 1], 'YTick', 0:0.1:1.0, ...
        'XTick', ticks, 'XTickLabel', sprintf('%d\n',subjectIndices(idx)), 'FontSize', 10);
    %set(gca,'TickLength',[0.002, 0.002])
    xtickangle(60)
    xlabel('subject id', 'FontSize', 14);
    ylabel('scores', 'FontSize', 14);
    grid on; box on;
    NicePlot.exportFigToPDF(sprintf('scoresAllMethods.pdf'), hFig, 300);
    fprintf('Done !\n');
    
    
    
    
    hFig = figure(12); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1830 370]);
    subplot('Position', [0.03 0.10 0.95 0.87]);
    
    plot(subjectIndices, rankingScores(idx),'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.5 0.8 1], 'MarkerSize', 10);
    hold on;
    plot(subjectIndices, otherScores(idx),  'ro-', 'LineWidth', 1.5, 'MarkerFaceColor', [1 0.5, 0.5], 'MarkerSize', 10);
    
    kk = find(idx==98)
    plot(kk*[1 1], [rankingScores(idx(kk)) otherScores(idx(kk))], 'ks-', 'MarkerSize', 18, 'LineWidth', 3);
    
    kk = find(idx==132);
    plot(kk*[1 1], [rankingScores(idx(kk)) otherScores(idx(kk))], 'ks-', 'MarkerSize', 18, 'LineWidth', 3);
    
    kk = find(idx==54);
    plot(kk*[1 1], [rankingScores(idx(kk)) otherScores(idx(kk))], 'ks-', 'MarkerSize', 18, 'LineWidth', 3);
    
    kk = find(idx==194);
    plot(kk*[1 1], [rankingScores(idx(kk)) otherScores(idx(kk))], 'ks-', 'MarkerSize', 18, 'LineWidth', 3);
    
    kk = find(idx==21);
    plot(kk*[1 1], [rankingScores(idx(kk)) otherScores(idx(kk))], 'ks-', 'MarkerSize', 18, 'LineWidth', 3);
    
    set(gca, 'XTick', ticks, 'XTickLabel', {}, 'YTick', 0:0.1:1, 'FontSize', 26);
    hL = legend({'PSF score', 'MTF score'});
    set(hL, 'FontSize', 18);
    xlabel('\it subject index', 'FontWeight', 'normal');
    ylabel('\it score', 'FontWeight', 'bold');
    NicePlot.exportFigToPDF(sprintf('ThibosRankingFigure.pdf'), hFig, 300);
end

% Method to compute the PSF matching scores
function plotSPectralWeights(wavelengths, psfSpectralWeights, otfSpectralWeights)

    hFig = figure(3); clf
    set(hFig, 'Position', [10 10 700 340], 'Color', [1 1 1]);
    subplot(1,2,1);
    stem(wavelengths, psfSpectralWeights, 'filled');
    set(gca, 'YLim', [-0.05 1.05], 'FontSize', 14);
    set(gca, 'XLim', [390 710]);
    axis 'square'; box on; grid on;
    xlabel('wavelength (nm)');
    ylabel('PSF weighting');
    subplot(1,2,2);
    stem(wavelengths, otfSpectralWeights, 'filled');
    set(gca, 'YLim', [-0.05 1.05], 'FontSize', 14);
    set(gca, 'XLim', [390 710]);
    axis 'square'; box on; grid on;
    xlabel('wavelength (nm)');
    ylabel('OTF weighting');
    drawnow     
end

function plotRankedSubjectData(subjectIndex, wavelengthsListToCompute, ...
    xSfCyclesDeg, xMinutes, meanZcoeffOTFSlicesX, meanSubjectOTFSlicesX, subjectOTFSlicesX, ...
    meanZcoeffOTFSlicesY, meanSubjectOTFSlicesY, subjectOTFSlicesY, ...
    rankedSubjectPSF, meanZcoeffPSF)

    nWaves = size(subjectOTFSlicesX,3);
    waveStride = 6;
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
    set(hFig, 'Position', [10 10 945 600], 'Color', [1 1 1]);
    
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