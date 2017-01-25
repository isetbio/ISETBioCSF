function varargout = v_IBIOCDColorGabor(varargin)
% varargout = v_IBIOCDColorGabor(varargin)
%
% Works by running t_colorGabor with various arguments and comparing
% results with those stored.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    %% Hello
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOCDColorGabor *****');
    
    %% Basic validation
    validationData1 = t_colorGabor('generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('vData1',validationData1, ...
        'UsingTheFollowingVariableTolerancePairs', ...
        'vData1.maxIsomerizations', 1e-5, ...
        'vData1.minIsomerizations', 1e-5);
    
end



