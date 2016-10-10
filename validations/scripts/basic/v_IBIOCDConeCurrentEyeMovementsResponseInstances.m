function varargout = v_IBIOCDConeCuurentEyeMovementsResponseInstances(varargin)
% varargout = v_IBIOCDConeCuurentEyeMovementsResponseInstances(varargin)
%
% Works by running t_coneCurrentEyeMovementsResponseInstances with various arguments and comparing
% results with those stored.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Hello
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOCDConeCuurentEyeMovementsResponseInstances *****');
    
    %% Basic validation
%     validationData1 = t_coneCurrentEyeMovementsResponseInstances('generatePlots',runTimeParams.generatePlots);
%     UnitTest.validationData('validationData1',validationData1);
    
    %% Spot version
    validationData2 = t_coneCurrentEyeMovementsResponseInstancesSpot('generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData2',validationData2);
end



