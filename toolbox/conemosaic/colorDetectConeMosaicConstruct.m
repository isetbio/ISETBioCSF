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
% 12/8/16  npc       Update it after linearized os model.


if (strcmp(mosaicParams.conePacking, 'hex'))
    resamplingFactor = 6;
    centerInMM = [0.0 0.0];                    % mosaic eccentricity in MM - this should obey mosaicParams.eccentricityDegs, but it does not do so yet
    spatiallyVaryingConeDensity = true;        % constant spatial density (at the mosaic's eccentricity)

    theMosaic = coneMosaicHex(resamplingFactor, spatiallyVaryingConeDensity, ...
                       'center', centerInMM*1e-3, ...
               'spatialDensity', [0 mosaicParams.LMSRatio]' ...
            );
else
    % Construct a cone mosaic with rectangular cone packing
    theMosaic = coneMosaic();
    theMosaic.spatialDensity = [0 mosaicParams.LMSRatio]';
end

% Set the outer segment model
if strcmp(mosaicParams.osModel, 'Linear')
    theMosaic.os = osLinear();
end

% Set mosaic field of view.  In principle this would be as large as the
% stimulus, but space and time considerations may lead to it being smaller.
if (isfield(mosaicParams, 'fieldOfViewDegs'))
    if (isa(theMosaic, 'coneMosaicHex'))
        theMosaic.setSizeToFOVForHexMosaic(mosaicParams.fieldOfViewDegs)
        theMosaic.visualizeGrid();
    else
        theMosaic.setSizeToFOV(mosaicParams.fieldOfViewDegs);
    end
end

% integration time
if (isfield(mosaicParams, 'integrationTimeInSeconds'))
    warningIntegrationTimeInSeconds = 25/1000;
    if (mosaicParams.integrationTimeInSeconds > warningIntegrationTimeInSeconds)
        fprintf(2, '\n\n=====================================================================\n');
        fprintf(2, 'Setting coneMosaic.integrationTime to %2.3f\n', mosaicParams.integrationTimeInSeconds);
        fprintf(2, 'Warning: Setting the coneMosaic.integrationTime > %2.3f is\n', warningIntegrationTimeInSeconds);
        fprintf(2, 'NOT recommended if you are interested in photocurrent computations.\n');
        fprintf(2, '=====================================================================\n');
    end
    theMosaic.integrationTime = mosaicParams.integrationTimeInSeconds;
end

% outer-segment time step
if (isfield(mosaicParams, 'osTimeStep'))
    if (mosaicParams.osTimeStepInSeconds > 1/1000)
        error('Cannot set os.timeStepto %2.4f. It must be less t<= to 1/1000, and preferably <= 0.5/1000', mosaicParams.osTimeStepInSeconds);
    end
    theMosaic.os.timeStep = mosaicParams.osTimeStepInSeconds;
end

% isomerization noise
if (isfield(mosaicParams, 'isomerizationNoise'))
    theMosaic.noiseFlag = mosaicParams.isomerizationNoise;
end

% Outer segment noise
if (isfield(mosaicParams, 'osNoise'))
    theMosaic.os.noiseFlag = mosaicParams.osNoise;
end