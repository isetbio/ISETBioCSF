function varargout = v_IBIOCDColorGabor(varargin)
% varargout = v_IBIOCDColorGabor(varargin)
%
% Works by running t_colorGabor with various arguments and comparing
% results with those stored.
%
% Validate applyKernel method of tfe parent class.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    
    %% Basic validation
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOColorDetectGaborScene *****');
    validationData1 = t_colorGabor([],'generatePlots',runTimeParams.generatePlots);
    UnitTest.validationData('validationData1',validationData1);
    
end



