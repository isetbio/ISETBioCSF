function varargout = v_IBIOCDIBanksEtAlReplicate(varargin)
% varargout = v_IBIOCDIBanksEtAlReplicate(varargin)
%
% Works by running t_coneIsomerizationsMovie with various arguments and comparing
% results with those stored.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    %% Hello
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOCDIBanksEtAlReplicate *****');
    
    %% Freeze reandom number generator seed
    rng('default');
    
    %% Basic validation
    validationData1 = c_BanksEtAlReplicate('nTrainingSamples',100,'cyclesPerDegree',10,'luminances',340,'pupilDiamMm',2,'generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData1',validationData1);
    
    validationData2 = c_BanksEtAlReplicate('nTrainingSamples',100,'cyclesPerDegree',10,'luminances',340,'pupilDiamMm',4,'generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData2',validationData2);
    
end



