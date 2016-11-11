function dirname = paramsToMosaicDirName(obj,mosaicParams)
% pdirname = paramsToMosaicDirName(obj,mosaicParams)
% 
% Generate a directory names that captures the mosaic parameters.

if (~strcmp(mosaicParams.type,'Mosaic')) && (~strcmp(mosaicParams.type,'Mosaic_v2'))
    error('Incorrect parameter type passed');
end

if (strcmp(mosaicParams.type,'Mosaic'))
    dirname = sprintf('LMS%0.2f_%0.2f_%0.2f_mfv%0.1f_ecc%0.1f_%s_time%0.1f_in%d_osn%d',...
        mosaicParams.LMSRatio(1),mosaicParams.LMSRatio(2),mosaicParams.LMSRatio(3), ...
        mosaicParams.fieldOfViewDegs, mosaicParams.eccentricityDegs, mosaicParams.conePacking, ...
        1000*mosaicParams.integrationTimeInSeconds,mosaicParams.isomerizationNoise,mosaicParams.osNoise);
else
    if (isscalar(mosaicParams.fieldOfViewDegs))
        mosaicFOVdegs = mosaicParams.fieldOfViewDegs*[1 1];
    else
        mosaicFOVdegs = [mosaicParams.fieldOfViewDegs(1) mosaicParams.fieldOfViewDegs(2)];
    end
    dirname = sprintf('[MOSAIC]_%sPacking_ecc%0.1f_LMSdensities%0.2f_%0.2f_%0.2f_FOVdegs%0.1fx%0.1f_intTime%0.1f_osTime%0.2f_pNoise%d_osNoise%d_eyeMvmnts%d',...
        mosaicParams.conePacking, mosaicParams.eccentricityDegs, ...
        mosaicParams.spatialLMSDensities(1),mosaicParams.spatialLMSDensities(2),mosaicParams.spatialLMSDensities(3), ...
        mosaicFOVdegs(1), mosaicFOVdegs(2),  ...
        1000*mosaicParams.integrationTimeSecs, 1000*mosaicParams.osTimeStepSecs, ...
        mosaicParams.photonNoise, mosaicParams.osNoise, ~mosaicParams.eyesDoNotMove);
end