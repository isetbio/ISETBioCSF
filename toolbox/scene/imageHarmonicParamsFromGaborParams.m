function imageHarmonicParams = imageHarmonicParamsFromGaborParams(gaborParams,colorModulationParams)
% imageHarmonicParams = imageHarmonicParamsFromGaborParams(gaborParams,colorModulationParams)
%
% The imageHarmonic function in isetbio takes a parameters structure that
% is a little different from how we'd like to specify scene parameters for
% this project. This function produces a structure to be passed to
% imageHarmonic based on the gaborParams and colorModulationParams structures
% we are using here.

% Set up base parameters
imageHarmonicParams = gaborParams;

% Pull contrast into the structure
imageHarmonicParams.contrast = colorModulationParams.contrast;

% Computed parameters.  These convert numbers to a form used by underlying
% routines.
cyclesPerImage = gaborParams.fieldOfViewDegs*gaborParams.cyclesPerDegree;
gaussianStdDegs = FWHMToStd(gaborParams.gaussianFWHMDegs);
gaussianStdImageFraction = gaussianStdDegs/gaborParams.fieldOfViewDegs;
imageHarmonicParams.freq = cyclesPerImage;

% Set GaborFlag.  Make it negative for half-cosine instead of Gaussian
% window.

imageHarmonicParams.GaborFlag = gaussianStdImageFraction;
switch gaborParams.windowType
    case 'halfcos'
        imageHarmonicParams.GaborFlag = -imageHarmonicParams.GaborFlag;
    case 'Gaussian'
    otherwise
        error('Unknown windowType passed');
end

end