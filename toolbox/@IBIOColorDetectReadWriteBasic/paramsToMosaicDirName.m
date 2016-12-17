function dirname = paramsToMosaicDirName(obj,mosaicParams) 
% pdirname = paramsToMosaicDirName(obj,mosaicParams)
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

if (strcmp(mosaicParams.conePacking, 'hex'))
    mosaicParams.coneSpacingMicrons = nan;
end
dirname = sprintf('[MOSAIC]_%sPacking_coneDiamMicrons%0.4f_coneSepMicrons%0.4f_ecc%0.1f_LMSdensities%0.2f_%0.2f_%0.2f_FOVdegs%0.1fx%0.1f_intTimeMilliSecs%0.0f_photonNoise%s_osModel%s_osTimeStepMilliSecs%0.2f_osNoise%s',...
    mosaicParams.conePacking, mosaicParams.innerSegmentDiamMicrons, mosaicParams.coneSpacingMicrons, mosaicParams.eccentricityDegs, ...
    mosaicParams.LMSRatio(1),mosaicParams.LMSRatio(2),mosaicParams.LMSRatio(3), ...
    mosaicFOVdegs(1), mosaicFOVdegs(2), ...
    1000*mosaicParams.integrationTimeInSeconds, sprintf('%s%s',upper(mosaicParams.isomerizationNoise(1)), mosaicParams.isomerizationNoise(2:end)), ...
    sprintf('%s%s',upper(mosaicParams.osModel(1)), mosaicParams.osModel(2:end)), ...
    1000*mosaicParams.osTimeStepInSeconds, ...
    sprintf('%s%s',upper(mosaicParams.osNoise(1)), mosaicParams.osNoise(2:end)));
end
