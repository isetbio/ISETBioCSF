function oiParams = oiParamsGenerate(varargin)
% oiParams = oiParamsGenerate(varargin)
%
% Properties related to computing the retinal image
%
%  fieldOfViewDegs - Field of view computed
%  offAxis - Boolean, compute falloff of intensity with position
%  blur - Boolean, incorporate optical blurring
%  lens - Boolean, incorporate filtering by lens
%
% Other parameters that are needed for the mosaic are
%  fieldOfViewDegs - Field of view computed
%  timeStepInSeconds - Time step to compute responses on
%  integrationTimeInSeconds - Cone integration time.  Generally the same as time step
% These get set outside this routine as they are typically matched up to
% other parameters set elsewhere.
%
% See also
%   responseParametersGenerate

oiParams.type = 'OpticsParams';

oiParams.offAxis = false;
oiParams.blur = true;
oiParams.lens = true;