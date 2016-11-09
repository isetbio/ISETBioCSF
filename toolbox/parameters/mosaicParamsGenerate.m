function mosaicParams = mosaicParamsGenerate(varargin)
% mosaicParams = mosaicParamsGenerate(varargin)
%
% Properties of the cone mosaic set here.
%  macular - Boolean, include macular pigment (may not be implemeted yet)
%  LMSRatio - Vector giving L, M, S cone ratio
%  isomerizationNoise - Boolean, add isomerization Poisson noise?
%  osModel - What model to use to compute photocurrent
%  conePacking - Type of cone mosaic
%    'RECT' - rectangular
%    'HEX' - hexagonal-like
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

mosaicParams.type = 'Mosaic';

mosaicParams.conePacking = 'RECT';
mosaicParams.macular = true;
mosaicParams.LMSRatio = [0.62 0.31 0.07];
mosaicParams.eccentricityDegs = 0;
mosaicParams.isomerizationNoise = false;
mosaicParams.osModel = 'Linear';
mosaicParams.osNoise = true;
