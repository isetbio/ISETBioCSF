function varargout = v_IBIOCDConeCuurentEyeMovementsMovie(varargin)
% varargout = v_IBIOCDConeCuurentEyeMovementsMovie(varargin)
%
% Works by running t_coneCuurentEyeMovementsMovie with various arguments and comparing
% results with those stored.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Hello
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOCDConeCuurentEyeMovementsMovie *****');
    
    %% Basic validation
    validationData1 = t_coneCurrentEyeMovementsMovie([],'generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData1',validationData1);
    
    %% Spot version
    validationData2 = t_coneCurrentEyeMovementsMovieSpot([],'generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData2',validationData2);
    
end



