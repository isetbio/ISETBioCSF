function varargout = v_IBIOCDColorDetectFindPerformance(varargin)
% varargout = v_IBIOCDColorDetectFindPerformance(varargin)
%
% Works by running t_colorDetectFindPerformance with various arguments and comparing
% results with those stored.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Hello
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOCDColorDetectFindPerformance *****');
    
    %% Basic validation
    validationData1 = t_colorDetectFindPerformance('generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData1',validationData1);
end



