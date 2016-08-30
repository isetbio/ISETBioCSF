function theBipolarMosaic = colorDetectBipolarMosaicConstruct(theMosaic, bipolarParams)
%  theBipolarMosaic = colorDetectBipolarMosaicConstruct(theMosaic, bipolarParams)
% 
% Construct bipolar mosaic, for use in IBIOColorDetect
% 
% See bipolarParamsGenerate for parameters including rectification type.
% 
% 8/2016 JRG (c) isetbio team

theBipolarMosaic = bipolar(theMosaic, bipolarParams);

theBipolarMosaic.set('sRFcenter',1);
theBipolarMosaic.set('sRFsurround',1);