function varargout = v_IBIOCDConeIsomerizationsMovie(varargin)
% varargout = v_IBIOCDConeIsomerizationsMovie(varargin)
%
% Works by running t_coneIsomerizationsMovie with various arguments and comparing
% results with those stored.
%
% Validate applyKernel method of tfe parent class.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    %% Freeze rng
    rng(1);
    
    %% Basic validation
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOColorDetectGaborScene *****');
    validationData1 = t_coneIsomerizationsMovie([],'generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData1',validationData1);
    
    %% Spot version
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOColorDetectGaborScene *****');
    validationData2 = t_coneIsomerizationsMovieSpot([],'generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData2',validationData2);
    
end



