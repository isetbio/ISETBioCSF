function [theOI, varargout] = colorDetectOpticalImageConstruct(oiParams, availableCustomWvfOpticsModels)
%colorDetectOpticalImageConstruct   Construct optical image object for use color detect simulations
%   theOI = colorDetectOpticalImageConstruct(oiParams, availableCustomWvfOpticsModels)
%
%   Construct optical image object, for use in the color detect simulations.
%
%   The passed parameters structure controls features of the optical
%   image object, allowing us to explore the effect of these features on performance.
%     oiParams.fieldOfViewDegrees - specify field of view.
%     oiParams.offAxis: true/false - do off axis vignetting or not
%     oiParams.blur: true/false - do optical blurring or not
%     oiParams.lens: true/false - put in human lens transmittance or not
%     oiParams.opticsModel - which model of human optics to use.

% 7/8/16  dhb  Wrote it.
% 6/2/17  npc  Added WvfHumanMeanOTFmagMeanOTFphase', 'WvfHumanMeanOTFmagZeroOTFphase', 'WvfHumanSingleSubject' optics.

varargout = {};

% Basic create and set field of view.
theOI = oiCreate('wvf human', oiParams.pupilDiamMm);

% Take opticsModel into account.
switch (oiParams.opticsModel)
    case 'WvfHuman'
    case availableCustomWvfOpticsModels
        fprintf('Computing custom OTF optics\n')
        [theOIdef, theCustomOI, Zcoeffs] = oiWithCustomOptics(oiParams.pupilDiamMm, oiParams.opticsModel);
        varargout{1} = Zcoeffs;
        plotOIs(theOIdef, theCustomOI);
        theOI = theCustomOI;
    case {'Geisler', 'GeislerLsfAsPsf', 'DavilaGeisler', 'DavilaGeislerLsfAsPsf', 'Westheimer', 'Williams'}
        theOI = ptb.oiSetPtbOptics(theOI,'opticsModel',oiParams.opticsModel);
    otherwise
        error('Unknown opticsModel string passed');
end

theOI = oiSet(theOI,'h fov',oiParams.fieldOfViewDegs);

% Set the pupil diamter
focalLength = oiGet(theOI,'distance');
desiredFNumber = focalLength/(oiParams.pupilDiamMm/1000);
theOI  = oiSet(theOI ,'optics fnumber',desiredFNumber);
pupilDiamMmCheck = 1000*oiGet(theOI,'optics aperture diameter');
if (max(abs(pupilDiamMmCheck - oiParams.pupilDiamMm)) > 1e-8)
    error('Failed to set pupil diameter as expected');
end

% Take out off axis vignetting if requested
optics = oiGet(theOI,'optics');
if (~oiParams.offAxis)
    optics = opticsSet(optics,'off axis method','skip');
end
theOI = oiSet(theOI,'optics',optics);

% Take out optical blurring if requested
optics = oiGet(theOI,'optics');
if (~oiParams.blur)
    optics = opticsSet(optics,'OTF',ones(size(opticsGet(optics,'OTF'))));
end
theOI = oiSet(theOI,'optics',optics);

% Take out lens transmittance if desired
if (~oiParams.lens)
    lens = oiGet(theOI,'lens');
    lens.density = 0;
	theOI = oiSet(theOI,'lens',lens);
end

end

function plotOIs(theOI, theCustomOI)
    figure(99); clf;
    
    wavelengthsList = [450:50:700];
    
    for k = 1:numel(wavelengthsList)
        wavelength = wavelengthsList(k);
        
        [otf, otf_fx, otf_fy, psf, psf_x, psf_y] = getOtfPsfData(theOI, wavelength);
        [theCustomOTF, otf_fx2, otf_fy2, theCustomPSF, psf_x2, psf_y2] = getOtfPsfData(theCustomOI, wavelength);

        if (any(otf_fx~=otf_fx2 | otf_fy~=otf_fy2  | psf_x~=psf_x2  | psf_y ~= psf_y2))
           error('Internal coordinate system transformation error');
        end

        %[k min(psf(:)) max(psf(:)) sum(psf(:)) min(theCustomPSF(:)) max(theCustomPSF(:)) sum(theCustomPSF(:))]

        subplot(2, numel(wavelengthsList), k)
        imagesc(psf_x, psf_y, psf);
        axis 'image';
        title(sprintf('%d nm', wavelength));
        subplot(2, numel(wavelengthsList), numel(wavelengthsList)+k)
        imagesc(psf_x, psf_y, theCustomPSF);
        axis 'image';
        colormap(gray(1024));
        drawnow;
    end
end

