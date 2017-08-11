function theMosaic = colorDetectConeMosaicConstruct(mosaicParams, varargin)
% colorDetectConeMosaicConstruct  Construct cone mosaic according to parameters structure
%   theMosaic = colorDetectConeMosaicConstruct(mosaicParams)
% 
% Construct a cone mosaic according to the passed parameters structure.
% Designed to allow us to control exactly what features of early vision
% we're using.
% 
%
% Key/Value pairs
%       'visualizeMosaic' - true/false (default true). Wether to visualize the cone mosaic
%
%
% 07/9/16  npc, dhb  Wrote it.
% 12/8/16  npc       Update it after linearized os model.

%% Parse arguments
p = inputParser;
p.addParameter('visualizeMosaic',true, @islogical);   
p.parse(varargin{:});
visualizeMosaic = p.Results.visualizeMosaic;

if (isfield(mosaicParams, 'fovDegs'))
    error('mosaicParams has fovDegs');
end

if (ischar(mosaicParams.conePacking))
    if (strcmp(mosaicParams.conePacking, 'hex'))
        resamplingFactor = 6;

        if (~isfield(mosaicParams, 'name'))
            mosaicParams.name = 'noname';
        end
        
        if (~isfield(mosaicParams, 'sConeMinDistanceFactor'))
            mosaicParams.sConeMinDistanceFactor = [];
        end
        
        if (~isfield(mosaicParams, 'sConeFreeRadiusMicrons'))
            mosaicParams.sConeFreeRadiusMicrons = [];
        end
        
        if (~isfield(mosaicParams, 'latticeAdjustmentPositionalToleranceF'))
            mosaicParams.latticeAdjustmentPositionalToleranceF = 0.01;
        end
        
        if (~isfield(mosaicParams, 'latticeAdjustmentDelaunayToleranceF'))
            mosaicParams.latticeAdjustmentDelaunayToleranceF = 0.001;
        end
        
        theMosaic = coneMosaicHex(resamplingFactor, ...
            'name', mosaicParams.name, ...
            'fovDegs', mosaicParams.fieldOfViewDegs, ...
            'eccBasedConeDensity', true, ...
            'sConeMinDistanceFactor', mosaicParams.sConeMinDistanceFactor, ...
            'sConeFreeRadiusMicrons', mosaicParams.sConeFreeRadiusMicrons, ...              
            'spatialDensity', [0 mosaicParams.LMSRatio]', ...
            'latticeAdjustmentPositionalToleranceF', mosaicParams.latticeAdjustmentPositionalToleranceF, ...   
            'latticeAdjustmentDelaunayToleranceF', mosaicParams.latticeAdjustmentDelaunayToleranceF ...
        );
        theMosaic.displayInfo();
        if (visualizeMosaic)
            theMosaic.visualizeGrid(); % 'visualizedConeAperture', 'geometricArea');
        end
    elseif (strcmp(mosaicParams.conePacking, 'hexReg'))
        resamplingFactor = 6;

        if (~isfield(mosaicParams, 'name'))
            mosaicParams.name = 'noname';
        end
        
        if (~isfield(mosaicParams, 'sConeMinDistanceFactor'))
            mosaicParams.sConeMinDistanceFactor = [];
        end
        
        if (~isfield(mosaicParams, 'sConeFreeRadiusMicrons'))
            mosaicParams.sConeFreeRadiusMicrons = [];
        end
        
        theMosaic = coneMosaicHex(resamplingFactor, ...
            'fovDegs', mosaicParams.fieldOfViewDegs, ...
            'customLambda', mosaicParams.coneSpacingMicrons, ...
            'customInnerSegmentDiameter', mosaicParams.innerSegmentSizeMicrons, ... 
            'eccBasedConeDensity', false, ...
            'sConeMinDistanceFactor', mosaicParams.sConeMinDistanceFactor, ...
            'sConeFreeRadiusMicrons', mosaicParams.sConeFreeRadiusMicrons, ...    
            'spatialDensity', [0 mosaicParams.LMSRatio]', ...
            'rotationDegs', mosaicParams.mosaicRotationDegs ...
        );
        theMosaic.displayInfo();
        
        if (visualizeMosaic)
            theMosaic.visualizeGrid(); % 'visualizedConeAperture', 'geometricArea');
        end
    elseif (strcmp(mosaicParams.conePacking, 'rect'))
        % Construct a cone mosaic with rectangular cone packing
        theMosaic = coneMosaic();
        theMosaic.spatialDensity = [0 mosaicParams.LMSRatio]';
        
        % Set the pigment collecting area   
        theMosaic.pigment.pdWidth = mosaicParams.innerSegmentSizeMicrons*1e-6;
        theMosaic.pigment.pdHeight = mosaicParams.innerSegmentSizeMicrons*1e-6;
        
        % Set the pigment geometric dimensions
        theMosaic.pigment.width = mosaicParams.coneSpacingMicrons*1e-6;
        theMosaic.pigment.height = mosaicParams.coneSpacingMicrons*1e-6;
        theMosaic.setSizeToFOV(mosaicParams.fieldOfViewDegs);
    else
        error('Unknown conePacking value:'' %s''.', mosaicParams.conePacking);
    end
else
    mosaicParams.conePacking
    error('Unknown conePacking value: ''%s''', mosaicParams.conePacking);
end

% Set whether to blur by cone aperture
theMosaic.apertureBlur = mosaicParams.apertureBlur;

% Set dark noise
theMosaic.coneDarkNoiseRate = mosaicParams.coneDarkNoiseRate;

% Set the outer segment model
if strcmp(mosaicParams.osModel, 'Linear')
    theMosaic.os = osLinear();
end

% Integration time
if (isfield(mosaicParams, 'integrationTimeInSeconds'))
    theMosaic.integrationTime = mosaicParams.integrationTimeInSeconds;
end

% Outer-segment time step
if (isfield(mosaicParams, 'osTimeStep'))
    if (mosaicParams.osTimeStepInSeconds > 1/1000)
        error('Cannot set os.timeStepto %2.4f. It must be less t<= to 1/1000, and preferably <= 0.5/1000', mosaicParams.osTimeStepInSeconds);
    end
    theMosaic.os.timeStep = mosaicParams.osTimeStepInSeconds;
end

% Isomerization noise
if (isfield(mosaicParams, 'isomerizationNoise'))
    theMosaic.noiseFlag = mosaicParams.isomerizationNoise;
end

% Outer segment noise
if (isfield(mosaicParams, 'osNoise'))
    theMosaic.os.noiseFlag = mosaicParams.osNoise;
end