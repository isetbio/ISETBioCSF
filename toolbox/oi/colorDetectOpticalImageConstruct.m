function theOI = colorDetectOpticalImageConstruct(oiParams)
%colorDetectOpticalImageConstruct   Construct optical image object for use color detect simulations
%   theOI = colorDetectOpticalImageConstruct(oiParams)
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

% Basic create and set field of view.
theOI = oiCreate('wvf human');

% Take opticsModel into account.
switch (oiParams.opticsModel)
    case 'WvfHuman'
    case 'WvfHumanMeanOTF'
        theCustomOI = updateOIWithCustomOptics(theOI, oiParams.pupilDiamMm, 'opticsModel',oiParams.opticsModel);
        plotOIs(theOI, theCustomOI);
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
fprintf('In plot OIs\n');
pause();
end

