function oiParams = oiParamsGenerate(varargin)
%oiParams  Properties related to computing the retinal image
% oiParams = oiParamsGenerate(varargin)
%
% This function stores parameters relevant to the optics model used in the
% calculations.
%
%  fieldOfViewDegs - Field of view computed
%  offAxis - Boolean, compute falloff of intensity with position
%  blur - Boolean, incorporate optical blurring
%  lens - Boolean, incorporate filtering by lens
%  opticsModel - String, what optics model to use
%    'WvfHuman'           Default.  Isetbio standard wavefront based model of human optics.
%    'DavilaGeisler'      PSF based on DavilaGeisler line spread function.
%    'Westheimer'         PSF based on Westheimer line spread function.
%    'Williams'           PSF based on Williams et al. MTF
%
% Other parameters that are needed for the mosaic are
%  fieldOfViewDegs - Field of view computed
%  integrationTimeInSeconds - Cone integration time.  Generally the same as time step
% These get set outside this routine as they are typically matched up to
% other parameters set elsewhere.
%
% See also
%   responseParametersGenerate, 

oiParams.type = 'Optics';

oiParams.offAxis = false;
oiParams.blur = true;
oiParams.lens = true;
oiParams.pupilDiamMm = 3;
oiParams.umPerDegree = 300;
oiParams.wavefrontSpatialSamples = 201;
oiParams.opticsModel = 'WvfHuman';
oiParams.opticalImagePadSizeDegs = [];