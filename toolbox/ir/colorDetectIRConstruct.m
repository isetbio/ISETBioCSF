function theIR = colorDetectIRConstruct(theBipolarMosaic, irParams)


% Create RGC object
theIR = ir(theBipolarMosaic, irParams);
% 
cellType = 'onParasol'; 
model = 'GLM';
theIR.mosaicCreate('type',cellType,'model',model);
% innerRetinaSU.mosaicCreate

% nTrials = 4; theIR = irSet(theIR,'numberTrials',nTrials);



