function imageHarmonicParams = imageHarmonicParamsFromGaborParams(gaborParams)
% imageHarmonicParams = imageHarmonicParamsFromGaborParams(gaborParams)
%
% The imageHarmonic function in isetbio takes a parameters structure that
% is a little different from how we'd like to specify scene parameters for
% this project. This function produces a structure to be passed to
% imageHarmonic based on the gaborParams structure we are using here.

% Set up base parameters
imageHarmonicParams = gaborParams;

% Computed parameters.  These convert numbers to a form used by underlying
% routines.
cyclesPerImage = gaborParams.fieldOfViewDegs*gaborParams.cyclesPerDegree;
gaussianStdDegs = FWHMToStd(gaborParams.gaussianFWHMDegs);
gaussianStdImageFraction = gaussianStdDegs/gaborParams.fieldOfViewDegs;
imageHarmonicParams.freq = cyclesPerImage;
imageHarmonicParams.GaborFlag = gaussianStdImageFraction;