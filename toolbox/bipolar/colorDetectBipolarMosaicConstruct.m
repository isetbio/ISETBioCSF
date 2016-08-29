function theBipolarMosaic = colorDetectBipolarMosaicConstruct(theMosaic, bipolarParams)

theBipolarMosaic = bipolar(theMosaic, bipolarParams);

theBipolarMosaic.set('sRFcenter',1);
theBipolarMosaic.set('sRFsurround',1);