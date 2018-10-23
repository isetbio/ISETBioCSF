function dirname = paramsToMosaicDirName(obj,mosaicParams) 
% dirname = paramsToMosaicDirName(obj,mosaicParams)
% 
% Generate a directory names that captures the mosaic parameters.

if (~strcmp(mosaicParams.type,'Mosaic'))
    error('Incorrect parameter type passed');
end
 
if (isscalar(mosaicParams.fieldOfViewDegs))
    mosaicFOVdegs = mosaicParams.fieldOfViewDegs*[1 1];
else
    mosaicFOVdegs = [mosaicParams.fieldOfViewDegs(1) mosaicParams.fieldOfViewDegs(2)];
end

% For hex mosaic, spacing is computed from retinal location, so 
% coneSpacingMicrons parameter has no effect.
if (strcmp(mosaicParams.conePacking, 'hex'))
    mosaicParams.coneSpacingMicrons = NaN;
end

if (mosaicParams.eccBasedConeQuantalEfficiency)
    mosaicParams.efficiencyCorrection = 'VariedConeEff';
else
    mosaicParams.efficiencyCorrection = 'FovealConeEff';
end

if (mosaicParams.eccBasedMacularPigment)
    mosaicParams.macularPigmentCorrection = 'VariedMacPigment';
else
    mosaicParams.macularPigmentCorrection = 'FovealMacPigment';
end

% For rect mosaic, rotation has no effect
if (strcmp(mosaicParams.conePacking, 'rect'))
    mosaicParams.mosaicRotationDegs = NaN;
end

% must keep this name short (for the file system)
%   255 chars for mac (hfs) / linux (ext)
%   242 for docker (aufs)
dirname = sprintf('M_%sPacking_coneSizeUm%0.4f_coneSepUm%0.4f_%s_%s_rotationDegs%d_eccentricityDegs%0.2f_LMSdensities%0.2f_%0.2f_%0.2f_FOVdegs%0.2fx%0.2f_intTimeMs%0.0f_photonNoise%s_osModel%s_osTimeStepMs%0.2f_osNoise%s_apBlur%d_dark_%d_%d_%d',...
    mosaicParams.conePacking, mosaicParams.innerSegmentSizeMicrons, mosaicParams.coneSpacingMicrons, ...
    mosaicParams.efficiencyCorrection, mosaicParams.macularPigmentCorrection, mosaicParams.mosaicRotationDegs,mosaicParams.eccentricityDegs, ...
    mosaicParams.LMSRatio(1),mosaicParams.LMSRatio(2),mosaicParams.LMSRatio(3), ...
    mosaicFOVdegs(1), mosaicFOVdegs(2), ...
    1000*mosaicParams.integrationTimeInSeconds, sprintf('%s%s',upper(mosaicParams.isomerizationNoise(1)), mosaicParams.isomerizationNoise(2:end)), ...
    sprintf('%s%s',upper(mosaicParams.osModel(1)), mosaicParams.osModel(2:end)), ...
    1000*mosaicParams.osTimeStepInSeconds, ...
    sprintf('%s%s',upper(mosaicParams.osNoise(1)), mosaicParams.osNoise(2:end)), ...
    mosaicParams.apertureBlur,...
    mosaicParams.coneDarkNoiseRate(1),...
    mosaicParams.coneDarkNoiseRate(2),...
    mosaicParams.coneDarkNoiseRate(3));
end
