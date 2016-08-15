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
% routines.  This one is frequency
cyclesPerImage = gaborParams.fieldOfViewDegs*gaborParams.cyclesPerDegree;
imageHarmonicParams.freq = cyclesPerImage;
% 
% % Set GaborFlag to specify window.  Different conventions about width for
% % half-cosine and Gaussian, plus make it negative for half-cosine instead of Gaussian
% % window.
switch gaborParams.windowType
    case 'halfcos'
        % Input specifies half the width of the half-cosine, and
        % imageHarmonic uses this convention for half-cosines.
        cosHalfWidthDegs = gaborParams.gaussianFWHMDegs;
        cosHalfWidthImageFraction = cosHalfWidthDegs/gaborParams.fieldOfViewDegs;
        imageHarmonicParams.GaborFlag = -cosHalfWidthImageFraction;
    case 'Gaussian'
        % Input specifies FWHM for Gaussian, and convert to standard
        % deviation which is what imageHarmonic wants.
        gaussianStdDegs = FWHMToStd(gaborParams.gaussianFWHMDegs);
        gaussianStdImageFraction = gaussianStdDegs/gaborParams.fieldOfViewDegs;
        imageHarmonicParams.GaborFlag = gaussianStdImageFraction;

    otherwise
        error('Unknown windowType passed');
end

end