function imageHarmonicParams = imageHarmonicParamsFromGaborParams(spatialParams, colorModulationParams)
% imageHarmonicParams = imageHarmonicParamsFromGaborParams(spatialParams,colorModulationParams)
%
% The imageHarmonic function in isetbio takes a parameters structure that
% is a little different from how we'd like to specify scene parameters for
% this project. This function produces a structure to be passed to
% imageHarmonic based on the spatialParams and colorModulationParams structures
% we are using here.

% Set up base parameters
imageHarmonicParams = spatialParams;

% Pull contrast into the structure
if (isstruct(colorModulationParams)) && (isfield(colorModulationParams, 'contrast'))
    imageHarmonicParams.contrast = colorModulationParams.contrast;
else
    imageHarmonicParams.contrast = colorModulationParams;
end
% Computed parameters.  These convert numbers to a form used by underlying
% routines.  This one is frequency
cyclesPerImage = spatialParams.fieldOfViewDegs*spatialParams.cyclesPerDegree;
imageHarmonicParams.freq = cyclesPerImage;
% 
% % Set GaborFlag to specify window.  Different conventions about width for
% % half-cosine and Gaussian, plus make it negative for half-cosine instead of Gaussian
% % window.
switch spatialParams.windowType
    case 'halfcos'
        % Input specifies half the width of the half-cosine, and
        % imageHarmonic uses this convention for half-cosines.
        cosHalfWidthDegs = spatialParams.gaussianFWHMDegs;
        cosHalfWidthImageFraction = cosHalfWidthDegs/spatialParams.fieldOfViewDegs;
        imageHarmonicParams.GaborFlag = -cosHalfWidthImageFraction;
    case 'Gaussian'
        % Input specifies FWHM for Gaussian, and convert to standard
        % deviation which is what imageHarmonic wants.
        gaussianStdDegs = FWHMToStd(spatialParams.gaussianFWHMDegs);
        gaussianStdImageFraction = gaussianStdDegs/spatialParams.fieldOfViewDegs;
        imageHarmonicParams.GaborFlag = gaussianStdImageFraction;

    otherwise
        error('Unknown windowType passed');
end

end