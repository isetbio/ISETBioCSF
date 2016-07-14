function theOI = colorDetectOpticalImageConstruct(oiParams)
% theOI = colorDetectOpticalImageConstruct(oiParams)
%
% Construct optical image object, for use in the color detect simulations.
%
% The passed parameters structure controls features of the optical
% image object, allowing us to explore the effect of these features on performance.
%
%   oiParams.fieldOfViewDegrees - specify field of view.
%   oiParams.offAxis: true/false - do off axis vignetting or not
%   oiParams.blur: true/false - do optical blurring or not
%   oiParams.lens: true/false - put in human lens transmittance or not
%
% 7/8/16  dhb  Wrote it.

theOI = oiCreate('human');
theOI = oiSet(theOI,'h fov',oiParams.fieldOfViewDegs);

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


