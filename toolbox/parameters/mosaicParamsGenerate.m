function mosaicParams = mosaicParamsGenerate(varargin)
%MOSAICPARAMSGENERATE  Set structure with properties of the cone mosaic
%    mosaicParams = MOSAICPARAMSGENERATE(varargin)
%
% Properties set here:
%  conePacking - Type of cone mosaic
%    'rect' - rectangular
%    'hexReg' - hexagonal, regular spacing
%    'hex'  - hexagonal, eccentricity-varying spacing %  macular - Boolean, include macular pigment (may not be implemeted yet)
%  LMSRatio - Vector giving L, M, S cone ratio
%  innerSegmetnSizeMicrons - Linear size of a square cone light-collecting area in microns value for coneMosaic.pigment.pdWidth
%  apertureBlur - Blur by the cone aperture?
%  coneSpacingMicrons - Cone spacing in microns
%  mosaicRotationDegs - Rotation of mosaic, in degrees.  Only applies to hex and hexReg.
%  macular - Model macular pigment (true/false)
%  eccentricityDegrees - Retinal eccentricity of modeled mosaic
%  integrationTimeInSeconds - Integration time for cone isomerizations
%  osTimeStepInSeconds - Time step for outer segment calculations
%  isomerizationNoise - Noise model for isomerizations
%  osModel - Model used for outer segment calculations
%  osNoise - Noise model for isomerizations
%  darkNoiseRate - LMS cone dark noise rate in iso/sec.
%
% Other parameters that are needed for the mosaic are (but not specified here) are:
%  fieldOfViewDegs - Field of view. This is computed based on stimulus.
%
% See also RESPONSEPARAMETERSGENERATE

mosaicParams.type = 'Mosaic';

mosaicParams.conePacking = 'rect';
mosaicParams.eccBasedConeQuantalEfficiency = false;
mosaicParams.LMSRatio = [0.62 0.31 0.07];
mosaicParams.innerSegmentSizeMicrons = 1.4;       
mosaicParams.apertureBlur = false;                  
mosaicParams.coneSpacingMicrons = 2.0;             
mosaicParams.mosaicRotationDegs = 0;
mosaicParams.macular = true;
mosaicParams.eccentricityDegs = 0;
mosaicParams.integrationTimeInSeconds = 10/1000;
mosaicParams.osTimeStepInSeconds = 0.1/1000;
mosaicParams.isomerizationNoise = 'none';           % Type coneMosaic.validNoiseFlags to get valid values
mosaicParams.osModel = 'Linear';
mosaicParams.osNoise = 'random';                    % Type outerSegment.validNoiseFlags to get valid values
mosaicParams.coneDarkNoiseRate = [0 0 0];

% Hex mosaic - specific params
mosaicParams.sConeMinDistanceFactor = 3.0;
mosaicParams.sConeFreeRadiusMicrons =  45;
mosaicParams.latticeAdjustmentPositionalToleranceF =  0.01;
mosaicParams.latticeAdjustmentDelaunayToleranceF = 0.001;
mosaicParams.maxGridAdjustmentIterations = Inf;
mosaicParams.marginF = [];
mosaicParams.resamplingFactor = 9;