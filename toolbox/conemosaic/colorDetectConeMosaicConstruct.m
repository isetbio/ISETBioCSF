function theMosaic = colorDetectConeMosaicConstruct(mosaicParams)
% theMosaic = colorDetectConeMosaicConstruct(mosaicParams)
% 
% Construct a cone mosaic according to the passed parameters structure.
% Designed to allow us to control exactly what features of early vision
% we're using.
% 
%   mosaicParams.fieldOfViewDegs - field of view in degrees
%   mosaicParams.LMSRatio - vector with three entries summing to one
%                           proportion of L, M, and S cones in mosaic
%
% THESE ARE NOT YET IMPLEMENTED
%   mosaicParams.macular -  true/false, include macular pigment?
%   mosaicParams.osModel - 'Linear','Biophys', which outer segment model
%
%  7/9/16  npc, dhb  Wrote it.
%

% Create a coneMosaic object here. 
if (isfield(mosaicParams, 'timeStepInSeconds'))
    if (isfield(mosaicParams, 'integrationTimeInSeconds'))
        theMosaic = coneMosaic(...
            'integrationTime',mosaicParams.integrationTimeInSeconds, ...  
                 'sampleTime',mosaicParams.timeStepInSeconds ...   % em, osCompute
                 );
    else
        theMosaic = coneMosaic(...
                 'sampleTime',mosaicParams.timeStepInSeconds ...   % em, osCompute
                 );
    end             
else
    if (isfield(mosaicParams, 'integrationTimeInSeconds'))
        theMosaic = coneMosaic('integrationTime',mosaicParams.integrationTimeInSeconds);
    else
        % Default foveal parameters for humans.
        theMosaic = coneMosaic;
    end
end


% Set mosaic field of view.  In principle this would be as large as the
% stimulus, but space and time considerations may lead to it being smaller.
if (isfield(mosaicParams, 'fieldOfViewDegs'))
    theMosaic.setSizeToFOV(mosaicParams.fieldOfViewDegs);
end

if (isfield(mosaicParams, 'photonNoise'))
    theMosaic.noiseFlag = mosaicParams.photonNoise;
end

if (isfield(mosaicParams, 'osNoise'))
    theMosaic.os.noiseFlag = mosaicParams.osNoise;
end


% Density of LMS cones
if (isfield(mosaicParams, 'LMSRatio'))
    if (numel(mosaicParams.LMSRatio) == 3)
        theMosaic.spatialDensity = [0 mosaicParams.LMSRatio(1) mosaicParams.LMSRatio(2) mosaicParams.LMSRatio(3)]';
    else
        theMosaic.spatialDensity = mosaicParams.LMSRatio(:);
    end
end