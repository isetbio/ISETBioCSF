function customOpticsPlayground()

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


function [theOI, theCustomOI, subjectOTFs] = generateCustomizedOI(wavelengthsToCompute, pupilDiameterMM, subjectsNum, magMethod, phaseMethods)

    % Load the 4.5 mm pupil diameter Thibos data set
    measPupilDiameterMM = 4.5;
    [Zcoeffs_SampleMean, Zcoeffs_S] = wvfLoadThibosVirtualEyes(measPupilDiameterMM);

    % Generate Z-coeffs for subjectsNum
    Zcoeffs = ieMvnrnd(Zcoeffs_SampleMean, Zcoeffs_S, subjectsNum)'; 
    
    [customMTFdataList, psfSingleSubject, subjectOTFs] = ...
        computeCustomMTF(Zcoeffs, Zcoeffs_SampleMean, pupilDiameterMM, wavelengthsToCompute, subjectsNum, magMethod, phaseMethods);
    
    theOI = oiCreate('wvf human', pupilDiameterMM);
    for k = 1:numel(phaseMethods)
        theCustomOI{k} = theOI;
        optics = oiGet(theCustomOI{k},'optics');
        wavelengthsListToCompute = opticsGet(optics,'wave');

        for waveIndex = 1:numel(wavelengthsListToCompute)
            theOtf = ifftshift(customMTFdataList{k}.otf);
            if (waveIndex == 1)
                otfData = zeros(size(theOtf,1), size(theOtf,2), numel(wavelengthsListToCompute));
            end
            otfData(:,:,waveIndex) = theOtf;
        end
        % Update optics with new OTF data
        optics = opticsSet(optics,'otf data',otfData);

        % Stick optics into oi
        theCustomOI{k} = oiSet(theCustomOI{k},'optics',optics);
    end
end

function [customMTFdataList, psfSingleSubject, otfSubjects] = computeCustomMTF(Zcoeffs, Zcoeffs_SampleMean, pupilDiameterMM, wavelengthsToCompute, subjectsNum, magMethod, phaseMethods)
    
    wvf_SampleMean = makeWVF(Zcoeffs_SampleMean, wavelengthsToCompute, pupilDiameterMM, 'sample mean');
    optics_SampleMean = oiGet(wvf2oi(wvf_SampleMean),'optics'); 
        
    otfSampleMean = opticsGet(optics_SampleMean, 'otf');
    otfSampleMeanMag = fftshift(abs(otfSampleMean));
    otfSampleMeanPhase = fftshift(angle(otfSampleMean));
    
    for subjectIndex = 1:subjectsNum
        % Make WVF for this subject
        wvfSubject = makeWVF(Zcoeffs(:,subjectIndex), wavelengthsToCompute, pupilDiameterMM, sprintf('subject-%d', subjectIndex));
        % Get PSF at selected wavelength
        psfSingleSubject = wvfGet(wvfSubject, 'psf', wavelengthsToCompute);
        optics = oiGet(wvf2oi(wvfSubject),'optics'); 
        
        % Get OTF at selected wavelength
        otf = opticsGet(optics, 'otf');
        sfx = opticsGet(optics, 'otf fx', 'cyclesperdeg');
        sfy = opticsGet(optics, 'otf fy', 'cyclesperdeg');
        otf = otfWithZeroCenteredPSF(otf, psfSingleSubject);
       
        if (subjectIndex == 1)
            otfSubjects = zeros(subjectsNum, size(otf,1), size(otf,2));
            customMTFdata.xSfGridCyclesDeg = sfx;
            customMTFdata.ySfGridCyclesDeg = sfy;
        end
        otfSubjects(subjectIndex,:,:) = otf;
    end % subjectIndex
    
    % Compute average MTF (mag and phase) across subjects
    otfMean = squeeze(mean(otfSubjects,1));
    otfMeanMag = fftshift(abs(otfMean));
    otfMeanPhase = fftshift(angle(otfMean));
    
    % Make radially symmetric magnitude
    if strcmp(magMethod, 'sample mean')
        magToUse = otfSampleMeanMag;
    elseif strcmp(magMethod, 'subject mean')
        magToUse = otfMeanMag;
    end
    
    centerPosition = floor(size(magToUse,1)/2) + 1;
    otfSlice = squeeze(magToUse(centerPosition, :));
    radiusMatrix = MakeRadiusMat(length(otfSlice),length(otfSlice),centerPosition);
    magToUse = interp1(0:floor(length(otfSlice)/2),otfSlice(centerPosition:end),radiusMatrix,'linear',0);
        
    for k = 1:numel(phaseMethods)
        % Compute MTF from radially symmetric magnitude and mean phase
        if strcmp(phaseMethods{k}, 'sample mean')
            phaseToUse = otfSampleMeanPhase;
        elseif strcmp(phaseMethods{k}, 'subject mean')
            phaseToUse = otfMeanPhase;
        elseif strcmp(phaseMethods{k}, 'zero')
            phaseToUse = 0*magToUse;
        end
        customMTFdata.otf = magToUse .* exp(1i*phaseToUse);
        customMTFdataList{k} = customMTFdata;
    end
end

function otf = otfWithZeroCenteredPSF(otf, psf)
    % Compute center of mass
    centerOfMass = computeCenterOfMass(psf);
   
    % Compute translation vector to bring center of mass at 0,0
    centerPosition = floor(size(psf,1)/2) + 1 * [1 1];
    translationVector = (centerPosition-centerOfMass);
    otf = shiftInFTplane(otf, translationVector);
end

function centerOfMass = computeCenterOfMass(psf0)
    [rc,cc] = ndgrid(1:size(psf0,1),1:size(psf0,2));
    Mt = sum(psf0(:));
    centerOfMassY = sum(psf0(:) .* rc(:)) / Mt;
    centerOfMassX = sum(psf0(:) .* cc(:)) / Mt;
    centerOfMass = [centerOfMassX centerOfMassY];
end

function theWVF = makeWVF(zcoeffs, wavelengthsToCompute, pupilDiameterMM, name)
    theWVF = wvfCreate(...
                'wave',wavelengthsToCompute,...
                'zcoeffs',zcoeffs,...
                'name', name);
    theWVF = wvfSet(theWVF,'calc pupil size',pupilDiameterMM);
    theWVF = wvfComputePSF(theWVF);
end

function [otf, otf_fx, otf_fy, psf, psf_x, psf_y] = getOtfPsfData(theOI, wavelengthsToCompute)
    theOptics = oiGet(theOI, 'optics');
    otf = abs(fftshift(opticsGet(theOptics,'otf data', wavelengthsToCompute)));
    otf_fx = opticsGet(theOptics, 'otf fx', 'cyclesperdeg');
    otf_fy = opticsGet(theOptics, 'otf fx', 'cyclesperdeg');

    psf = opticsGet(theOptics,'psf data', wavelengthsToCompute);
    %psfSupport = opticsGet(theOptics,'psf support', 'deg');
    
    [X,Y] = meshgrid(otf_fx,otf_fy);
    [xGridMinutes,yGridMinutes] = SfGridCyclesDegToPositionGridMinutes(X,Y);
    
    % Make it minutes of arc
    psf_x = xGridMinutes(1,:);
    psf_y = psf_x;
end


function otf = shiftInFTplane(otf, translationVector)
  
    [N, M] = size(otf);
    x = [0:floor(N/2) floor(-N/2)+1:-1];
    y = [0:floor(M/2) floor(-M/2)+1:-1];
    x_shift = exp(-1i * 2 * pi * translationVector(2) * x' / N);
    y_shift = exp(-1i * 2 * pi * translationVector(1) * y / M);

    % Force conjugate symmetry. Otherwise this frequency component has no
    % corresponding negative frequency to cancel out its imaginary part.
    if mod(N, 2) == 0
        x_shift(N/2+1) = real(x_shift(N/2+1));
    end 
    if mod(M, 2) == 0
        y_shift(M/2+1) = real(y_shift(M/2+1));
    end
    otf = otf .* (x_shift * y_shift);
end
