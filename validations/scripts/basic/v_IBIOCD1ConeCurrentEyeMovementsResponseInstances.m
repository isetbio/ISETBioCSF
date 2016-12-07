function varargout = v_IBIOCD1ConeCurrentEyeMovementsResponseInstances(varargin)
% varargout = v_IBIOCD1ConeCuurentEyeMovementsResponseInstances(varargin)
%
% Works by running t_coneCurrentEyeMovementsResponseInstances with various arguments and comparing
% results with those stored.
%
% The 1 in the filename is to make sure that's gets run in the right order
% on a run through all validation scripts.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Hello
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOCD1ConeCuurentEyeMovementsResponseInstances *****');
    
    %% Basic validation
    [validationData1, extraData1] = t_coneCurrentEyeMovementsResponseInstances('generatePlots',runTimeParams.generatePlots, 'freezeNoise', true);
    UnitTest.validationData('validationData1',validationData1);
    UnitTest.extraData('extraData1',extraData1);
     
    %% Spot version
    [validationData2, extraData2] = t_coneCurrentEyeMovementsResponseInstancesSpot('generatePlots',runTimeParams.generatePlots, 'freezeNoise', true);
    UnitTest.validationData('validationData2',validationData2);
    UnitTest.extraData('extraData2',extraData2);
end




