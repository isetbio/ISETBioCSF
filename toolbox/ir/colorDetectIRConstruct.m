function theIR = colorDetectIRConstruct(theBipolarMosaic, irParams)
% theIR = colorDetectIRConstruct(theBipolarMosaic, irParams)
% 
% Construct the inner retina object for IBIOColorDetect tutorials.
% 
% See irParamsGenerate for the parameters that can be altered, mainly
% eccentricity and polar angle of simulated retinal patch.
% 
% The cell type may be altered here.
% 
% 8/2016 JRG (c) isetbio team

% Create RGC object
theIR = ir(theBipolarMosaic, irParams);
% 
cellType = 'onParasol'; 
model = 'GLM';
theIR.mosaicCreate('type',cellType,'model',model);
% innerRetinaSU.mosaicCreate

% nTrials = 4; theIR = irSet(theIR,'numberTrials',nTrials);



