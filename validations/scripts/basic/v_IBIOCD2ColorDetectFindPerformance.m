function varargout = v_IBIOCD2ColorDetectFindPerformance(varargin)
% varargout = v_IBIOCDC2olorDetectFindPerformance(varargin)
%
% Works by running t_colorDetectFindPerformance with various arguments and comparing
% results with those stored.
%
% The 2 in the filename is to make sure that's gets run in the right order
% on a run through all validation scripts.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Hello
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOCD2ColorDetectFindPerformance *****');
    
    %% Basic validation
    [validationData1,extraData1] = t_colorDetectFindPerformance('generatePlots',runTimeParams.generatePlots,'plotPsychometric',true,'freezeNoise',true, 'employStandardHostComputerResources', true);
    UnitTest.validationData('validationData1',validationData1);
    UnitTest.extraData('extraData1',extraData1);
end



