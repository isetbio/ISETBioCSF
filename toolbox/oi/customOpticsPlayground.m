% Demonstrates mean Zcoeff vs meanOTF OTF and PSFs
function customOpticsPlayground()

    oiParams.fieldOfViewDegs = 2.0;
    oiParams.offAxis = false;

    pupilDiamMM = [2.0 3.0];
    wavelengthsToDisplay = [450 550 650];
    
    figNo = 0;
    for k = 1:numel(pupilDiamMM)
        oiParams.pupilDiamMM = pupilDiamMM(k);
        fprintf('Computing optics for %2.2fmm pupil diameter\n', oiParams.pupilDiamMM);
        % Compute mean Z-coeff based and mean OTF based optics
        [theMeanZcoeffOI, theMeanOTFBasedOI, Zcoeffs, populationOTFdata] = oiWithCustomOptics(oiParams.pupilDiamMM, 'WvfHumanMeanOTFmagMeanOTFphase');
        % Adjust OIs according to oiParams
        theMeanZcoeffOI = finalizeOI(theMeanZcoeffOI, oiParams);
        theMeanOTFBasedOI = finalizeOI(theMeanOTFBasedOI, oiParams);
        % Plot OTF/PSF for different wavelengths
        for l = 1:numel(wavelengthsToDisplay)
            figNo = figNo + 1;
            plotOI(figNo, 1, theMeanZcoeffOI,   populationOTFdata, wavelengthsToDisplay(l), 'mean Zcoeff', oiParams.pupilDiamMM);
            plotOI(figNo, 2, theMeanOTFBasedOI, populationOTFdata, wavelengthsToDisplay(l), 'mean OTF',    oiParams.pupilDiamMM);
        end % l
    end % k
end

function plotOI(figNo, figRow, theOI, populationOTFdata, wavelengthToDisplay, theOIlabel, pupilDiam)

    % Extract data for plotting
    [otf, otf_fx, otf_fy, psf, psf_x, psf_y] = getOtfPsfData(theOI, wavelengthToDisplay);
    centerRow = find(otf_fy == 0);
    otfSliceX = squeeze(otf(centerRow,:));
    wavelengthSupport = opticsGet(oiGet(theOI, 'optics'), 'wave');
    wIndex = find(wavelengthSupport == wavelengthToDisplay);
    populationOTFdata = abs(squeeze(populationOTFdata(:,:,:,wIndex)));
    
    % Set up figure
    hFig = figure(figNo); 
    if (figRow == 1)
        clf;
        set(hFig, 'Position', [10 10 1160 770], 'Color', [1 1 1]);
    end
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 3, ...
           'heightMargin',   0.1, ...
           'widthMargin',    0.1, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.04);
       
    % Plot the 2D OTF
    subplot('Position', subplotPosVectors(figRow,1).v);
    imagesc(otf_fx, otf_fy, otf);
    hold on;
    plot([otf_fx(1) otf_fx(end)], [0 0], 'r-');
    plot([0 0], [otf_fy(1) otf_fy(end)], 'r-');
    hold off;
    axis 'image';
    set(gca, 'XLim', [-60 60], 'YLim', [-60 60]);
    xlabel('spatial frequency (c/deg)');
    title(sprintf('pupilDiam: %2.1fmm, lambda: %d nm\n2D OTF (method: %s)', pupilDiam, wavelengthToDisplay, theOIlabel));
    set(gca, 'FontSize', 14);
    
    % Plot the OTF slices
    subplot('Position', subplotPosVectors(figRow,2).v);
    hold on
    for subjectIndex = 1:size(populationOTFdata, 1)
        subjectOTF = fftshift(squeeze(populationOTFdata(subjectIndex,:,:)));
        subjectSliceX  = squeeze(subjectOTF(centerRow,:));
        plot(otf_fx, subjectSliceX, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
    end
    plot(otf_fx, otfSliceX, 'g-', 'LineWidth', 1.5);
    hold off;
    axis 'square';
    set(gca, 'XLim', [0 60], 'YLim', [0 1], 'XTick', 0:10:100, 'YTick', 0:0.2:1.0);
    xlabel('spatial frequency (c/deg)');
    grid on
    title(sprintf('1D OTF (method: %s)', theOIlabel));
    set(gca, 'FontSize', 14);
    
    % Plot the PSF
    subplot('Position', subplotPosVectors(figRow,3).v);
    imagesc(psf_x, psf_y, psf);
    hold on;
    plot([psf_x(1) psf_x(end)], [0 0], 'r-');
    plot([0 0], [psf_y(1) psf_y(end)], 'r-');
    hold off;
    axis 'image';
    set(gca, 'XLim', 6*[-1 1], 'YLim', 6*[-1 1]);
    xlabel('x-space (arc min)');
    title(sprintf('2D PSF (method: %s)', theOIlabel));
    set(gca, 'FontSize', 14);
    
    colormap(gray(1024))
    drawnow;
end

    
function theOI = finalizeOI(theOI, oiParams)
    % Set the FOV
    theOI = oiSet(theOI,'h fov',oiParams.fieldOfViewDegs);
    
    % Set the f-number based on the pupil diameter
    focalLength = oiGet(theOI,'distance');
    desiredFNumber = focalLength/(oiParams.pupilDiamMM/1000);
    theOI  = oiSet(theOI ,'optics fnumber',desiredFNumber);
    
    % Adjust the off-axis fall-off
    optics = oiGet(theOI,'optics');
    if (~oiParams.offAxis)
        optics = opticsSet(optics,'off axis method','skip');
    end
    theOI = oiSet(theOI,'optics',optics);
end



function customOpticsPlaygroundOLD()

    phaseMethods = {'zero', 'subject mean'};  % 'subject mean' or 'zero'
    
    wavelengthsListToCompute = 450;
    pupilDiameterMM = 2.0;
    
    for k = 1:numel(wavelengthsListToCompute)
        doIt(phaseMethods, wavelengthsListToCompute(k), pupilDiameterMM);
    end
    
end

function doIt(phaseMethods, wavelengthsToCompute, pupilDiameterMM)
    
    subjectsNum = 200;
    magMethod = 'subject mean';
    [theOI, theCustomOIs, subjectOTFs] = generateCustomizedOI(wavelengthsToCompute, pupilDiameterMM, subjectsNum, magMethod, phaseMethods);
    [theOTF, otf_fx, otf_fy, thePSF, psf_x, psf_y] = getOtfPsfData(theOI, wavelengthsToCompute);
    [theCustomOTF{1}, otf_fx, otf_fy, theCustomPSF{1}, psf_x, psf_y] = getOtfPsfData(theCustomOIs{1}, wavelengthsToCompute);
    [theCustomOTF{2}, otf_fx, otf_fy, theCustomPSF{2}, psf_x, psf_y] = getOtfPsfData(theCustomOIs{2}, wavelengthsToCompute);
    
    theOptics = oiGet(theOI, 'optics');
    wls = opticsGet(theOptics,'wave');
    wIndex = find(wls == wavelengthsToCompute);
    
    selectedRow = floor(size(theOTF,1)/2) + 1;
    otfSliceRow = squeeze(theOTF(selectedRow,:));
    otfCustomSliceRow{1} = squeeze(theCustomOTF{1}(selectedRow,:));
    otfCustomSliceRow{2} = squeeze(theCustomOTF{2}(selectedRow,:));
    
    sfRange = [0 70];
    psfRange = [0 5];
    psfTicks = -5:1:5;
    if (wavelengthsToCompute >= 550) && (wavelengthsToCompute < 650)
        psfRange = [0 0.8];
        psfTicks = [-0.8:0.2:0.8];
    end
    if (wavelengthsToCompute < 450)
        psfRange = [0 3];
    end
    
    hFig = figure(10); clf;
    set(hFig, 'Position', [10 10 830 950], 'Color', [1 1 1]);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 1+numel(phaseMethods), ...
           'colsNum', 3, ...
           'heightMargin',   0.08, ...
           'widthMargin',    0.04, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.001, ...
           'bottomMargin',   0.05, ...
           'topMargin',      0.04);
    subplot('Position', subplotPosVectors(1,1).v);
    hold on
    for k = 1:size(subjectOTFs,1)
        otfSubjectSlice = abs(fftshift(squeeze(subjectOTFs(k,:,:))));
        otfSubjectSliceRow = squeeze(otfSubjectSlice(selectedRow,:));
        plot(otf_fx, otfSubjectSliceRow, 'k-', 'LineWidth', 1.0, 'Color', [0.5 0.5 0.5]);
    end
    plot(otf_fx, otfSliceRow, 'b-', 'LineWidth', 2.0);
    set(gca, 'XLim', sfRange, 'YLim', [0 1], 'FontSize', 12);
    title(sprintf('OTF slices (%d nm) \nisetbio default', wavelengthsToCompute));
    
    subplot('Position', subplotPosVectors(1,2).v);
    imagesc(otf_fx, otf_fy, theOTF);
    hold on;
    plot([0 otf_fx(end)], otf_fy(selectedRow)*[1 1], 'b-', 'LineWidth', 1.5);
    hold off
    axis 'image';
    set(gca, 'XLim', sfRange(2)*[-1 1], 'YLim', sfRange(2)*[-1 1], 'FontSize', 12);
    title(sprintf('OTF (%d nm)\nisetbio default', wavelengthsToCompute));
    
    subplot('Position', subplotPosVectors(1,3).v);
    imagesc(psf_x, psf_y, thePSF);
    hold on;
    psfMiddleRow = find(psf_y==0);
    psfSlice = squeeze(thePSF(psfMiddleRow, :));
    plot(psf_x, -psfRange(2) + 2*psfRange(2)*psfSlice / max(psfSlice), 'c-', 'LineWidth', 1.5);
    hold off;
    axis 'image'; axis 'xy'
    set(gca, 'XTick', psfTicks, 'YTick', psfTicks, 'XLim', psfRange(2)*[-1 1], 'YLim', psfRange(2)*[-1 1], 'FontSize', 12);
    title(sprintf('PSF (%d nm, sum: %2.4f)\nisetbio default', wavelengthsToCompute, sum(thePSF(:))));

    for phaseMethodIndex = 1: numel(phaseMethods)
    subplot('Position', subplotPosVectors(1+phaseMethodIndex,1).v);
    hold on
    for k = 1:size(subjectOTFs,1)
        otfSlice = abs(fftshift(squeeze(subjectOTFs(k,:,:))));
        otfSubjectSliceRow = squeeze(otfSlice(selectedRow,:));
        plot(otf_fx, otfSubjectSliceRow, 'k-', 'LineWidth', 1.0, 'Color', [0.5 0.5 0.5]);
    end
    plot(otf_fx, otfCustomSliceRow{phaseMethodIndex}, 'r-', 'LineWidth', 2.0);
    set(gca, 'XLim', sfRange, 'YLim', [0 1], 'FontSize', 12);
    if (phaseMethodIndex == numel(phaseMethods))
        xlabel('spatial frequency (c/deg)');
    end
    
    title(sprintf('OTF slices'));
    
    subplot('Position', subplotPosVectors(1+phaseMethodIndex,2).v);
    imagesc(otf_fx, otf_fy, theCustomOTF{phaseMethodIndex});
    hold on
    plot([0 otf_fx(end)], otf_fy(selectedRow)*[1 1], 'r-', 'LineWidth', 1.5);
    hold off
    axis 'image';
    if (phaseMethodIndex == numel(phaseMethods))
        xlabel('spatial frequency (c/deg)');
    end
    set(gca, 'XLim', sfRange(2)*[-1 1], 'YLim', sfRange(2)*[-1 1], 'FontSize', 12);
    title(sprintf('OTF \nmag:''%s'', phase: ''%s''', magMethod, phaseMethods{phaseMethodIndex}));
    
    subplot('Position', subplotPosVectors(1+phaseMethodIndex,3).v);
    psf = theCustomPSF{phaseMethodIndex};
    imagesc(psf_x, psf_y, psf);
    hold on;
    psfCustomSlice = squeeze(theCustomPSF{phaseMethodIndex}(psfMiddleRow, :));
    plot(psf_x, -psfRange(2) + 2*psfRange(2)*psfCustomSlice / max(psfCustomSlice), 'm-', 'LineWidth', 1.5);
    plot(psf_x, -psfRange(2) + 2*psfRange(2)*psfSlice / max(psfSlice), 'c-', 'LineWidth', 1.5);
    hold off;
    axis 'image'; axis 'xy'
    set(gca, 'XTick', psfTicks, 'YTick', psfTicks, 'XLim', psfRange(2)*[-1 1], 'YLim', psfRange(2)*[-1 1], 'FontSize', 12);
    title(sprintf('PSF (sum: %2.4f)', sum(psf(:))));
    end % phaseMethodIndex
    
    colormap(gray);
     
    drawnow
    NicePlot.exportFigToPDF(sprintf('Mag_%s_Wave_%dnm.pdf', magMethod, wavelengthsToCompute), hFig, 300);
end

