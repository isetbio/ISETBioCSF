function mosaicParams = mosaicParamsGenerate(varargin)
%mosaicParams  Set structure with properties of the cone mosaic
%    mosaicParams = mosaicParamsGenerate(varargin)
%
% Properties set here:
%  conePacking - Type of cone mosaic
%    'rect' - rectangular
%    'hexReg' - hexagonal, regular spacing
%    'hex'  - hexagonal, eccentricity-varying spacing %  macular - Boolean, include macular pigment (may not be implemeted yet)
%  LMSRatio - Vector giving L, M, S cone ratio
%  innerSegmetnSizeMicrons - Linear size of a square cone light-collecting area in microns value for coneMosaic.pigment.pdWidth
%  apertureBlur - blur by the cone aperture?
%  coneSpacingMicrons - Cone spacing in microns
%  macular - Model macular pigment (true/false)
%  eccentricityDegrees - Retinal eccentricity of modeled mosaic
%  integrationTimeInSeconds - Integration time for cone isomerizations
%  osTimeStepInSeconds - Time step for outer segment calculations
%  isomerizationNoise - Noise model for isomerizations
%  osModel - Model used for outer segment calculations
%  osNoise - Noise model for isomerizations
%
% Other parameters that are needed for the mosaic are (but not specified
% here) are:
%  fieldOfViewDegs - Field of view computed
%
% See also
%   responseParametersGenerate

mosaicParams.type = 'Mosaic';

mosaicParams.conePacking = 'rect';
mosaicParams.LMSRatio = [0.62 0.31 0.07];
mosaicParams.innerSegmentSizeMicrons = 1.4;       
mosaicParams.apertureBlur = false;                  % Blur by cone aperture?
mosaicParams.coneSpacingMicrons = 2.0;              % Cone spacing in microns, isetbio default value for coneMosaic.pigment.width
mosaicParams.macular = true;
mosaicParams.eccentricityDegs = 0;
mosaicParams.integrationTimeInSeconds = 10/1000;
mosaicParams.osTimeStepInSeconds = 0.1/1000;
mosaicParams.isomerizationNoise = 'none';           % Type coneMosaic.validNoiseFlags to get valid values
mosaicParams.osModel = 'Linear';
mosaicParams.osNoise = 'random';                    % Type outerSegment.validNoiseFlags to get valid values               