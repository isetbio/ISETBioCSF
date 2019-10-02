function generateMosaicWithPSFsuperimposedFig
% A function to display the mosaic with the superimposed PSF
%
% Syntax:
%   generateMosaidWithPSFsuperimposedFig
%
% Description:
%    Display the figure containing the mosaic with the superimposed PSF
%    contained in Cottaris et al 2019(A).
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    [rootPath, ~] = fileparts(which(mfilename));
    rootPath = strrep(rootPath, 'script', 'isetbio_resources');

    mosaicFOV = 0.6;
    load(fullfile(rootPath, ...
        sprintf('coneMosaic_%1.2fdegFOV.mat', mosaicFOV)), ...
        'theConeMosaic');

    hFig = figure(4);
    clf
    set(hFig, 'Position', [10 10 760 350], 'Color', [1 1 1]);
    ax = subplot(1, 2, 1);
    
    psfRangeArcMin = 4;
    
    theConeMosaic.visualizeGrid('axesHandle', ax, ...
        'visualizedConeAperture', 'geometricArea', ...
        'apertureShape', 'disks', 'backgroundColor', [1 1 1], ...
        'overlayHexMesh', true);

    xTicksArcMin = (-psfRangeArcMin:1:psfRangeArcMin)
    xTicksMeters = xTicksArcMin / 60 * ...
        theConeMosaic.micronsPerDegree * 1e-6;
    set(ax, 'XTickLabels', xTicksArcMin, 'XTick', xTicksMeters, ...
        'YTickLabels', xTicksArcMin, 'YTick', xTicksMeters);
    xlabel(ax, 'arc min');
    set(ax, 'XLim', [xTicksMeters(1) xTicksMeters(end)], ...
        'YLim', [xTicksMeters(1) xTicksMeters(end)]);
    set(gca, 'FontSize', 16);
    xlabel('arc min');
    ylabel('arc min');
    box(ax, 'on')
    sceneFOV = 1.5;
    opticsModel = 'ThibosAverageSubject3MMPupil';
    %opticsModel = 'ThibosBestPSFSubject3MMPupil';

    pupilDiamMm = 3.0;
    theOI = generateOI(opticsModel, pupilDiamMm, sceneFOV);

    micronsPerDegree = 300;
    optics = oiGet(theOI, 'optics');
    waves = opticsGet(optics, 'wave');
    psfSupportMicrons = opticsGet(optics, 'psf support', 'um');

    xGridMicrons = psfSupportMicrons{1};
    yGridMicrons = psfSupportMicrons{2};
    xGridMinutes = 60 * psfSupportMicrons{1} / micronsPerDegree;
    yGridMinutes = 60 * psfSupportMicrons{2} / micronsPerDegree;

    xx = find(abs(xGridMinutes) <= psfRangeArcMin);
    yy = find(abs(yGridMinutes) <= psfRangeArcMin);

    targetWavelength = 450;
    [~, idx] = min(abs(targetWavelength - waves));
    targetWavelength = waves(idx);

    psf = opticsGet(optics, 'psf data', targetWavelength);

    ax = subplot(1, 2, 2);
    contourLevels = 0:0.1:1.0;
    contourf(ax, xGridMinutes, yGridMinutes, ...
        psf / max(psf(:)), contourLevels);

    axis 'image';
    axis 'xy';
    set(gca, 'XLim', [xTicksArcMin(1) xTicksArcMin(end)], ...
        'YLim', [xTicksArcMin(1) xTicksArcMin(end)], 'CLim', [0 1], ...
            'XTick', xTicksArcMin, 'YTick', xTicksArcMin);
    set(gca, 'FontSize', 16);
    xlabel('arc min');
    title(sprintf('%s\n @%d nm', opticsModel, targetWavelength));
    colormap(1 - gray)

    NicePlot.exportFigToPDF(...
        sprintf('%s-%d nm.pdf', opticsModel, targetWavelength), hFig, 300)
end

function theOI = generateOI(opticsModel, pupilDiamMm, horizontalFOV)
% Generate an optical image
%
% Syntax:
%   theOI = generateOI(opticsModel, pupilDiamMm, horizontalFOV)
%
% Description:
%    Generate the optical image necessary for the figure.
%
% Inputs:
%    opticsModel   - String. The name of the optics model to use.
%    pupilDiamMm   - Numeric. The pupil diameter in millimeters.
%    horizontalFOV - Numeric. The horizontal field of view.
%
% Outputs:
%    theOI         - 
%
% Optional key/value pairs:
%    None.
%

    oiParams = struct(...
        'opticsModel', opticsModel, ...
        'wavefrontSpatialSamples', 301, ...
        'pupilDiamMm', pupilDiamMm, ...
        'umPerDegree', 300);
    theOI = oiWithCustomOptics(oiParams.opticsModel, ...
        oiParams.wavefrontSpatialSamples, ...
        oiParams.pupilDiamMm, oiParams.umPerDegree);

    % Set the FOV
    theOI = oiSet(theOI, 'h fov', horizontalFOV);

    % Set the fNumber
    focalLength = oiGet(theOI, 'distance');
    desiredFNumber = focalLength / (oiParams.pupilDiamMm / 1000);
    theOI  = oiSet(theOI , 'optics fnumber', desiredFNumber);
end
