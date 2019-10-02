function varargout = v_IBIOCDColorGabor(varargin)
% Validate the ISETBioCSF color gabor
%
% Syntax:
%   [varargout] = v_IBIOCDColorGabor([varargin])
%
% Description:
%    Works by running t_colorGabor with various arguments and comparing
%    results with those stored.
%
% Inputs:
%    None required.
%
% Outputs:
%    None required.
%
% Optional key/value pairs:
%    NEEDS TO BE ADDED.
%

    varargout = ...
        UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);

end

function ValidationFunction(runTimeParams)
% Implement the ISETBio validation code
%
% Syntax:
%   ValidationFunction(runTimeParams)
%
% Description:
%    Implement the ISETBio validation code.
%
% Inputs:
%    runTimeParams - Struct. A structure containing the runtime parameters.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
    %% Hello
    UnitTest.validationRecord('SIMPLE_MESSAGE', ...
        '***** v_IBIOCDColorGabor *****');

    %% Basic validation
    validationData1 = ...
        t_colorGabor('generatePlots', runTimeParams.generatePlots);
    UnitTest.validationData('vData1', validationData1, ...
        'UsingTheFollowingVariableTolerancePairs', ...
        'vData1.maxIsomerizations', 1e-5, ...
        'vData1.minIsomerizations', 1e-5);

end
