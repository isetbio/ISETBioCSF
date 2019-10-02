function varargout = v_IBIOCDConeIsomerizationsMovie(varargin)
% Validate the cone isomerizations movie
%
% Syntax:
%   [varargout] = v_IBIOCDConeIsomerizationsMovie([varargin])
%
% Description:
%    Works by running t_coneIsomerizationsMovie with various arguments and
%    comparing results with those stored.
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
        '***** v_IBIOCDConeIsomerizationsMovie *****');

    %% Basic validation
    validationData1 = t_coneIsomerizationsMovie('generatePlots', ...
        runTimeParams.generatePlots);
    UnitTest.validationData('validationData1', validationData1);

    %% Spot version
    validationData2 = t_coneIsomerizationsMovieSpot(...
        'generatePlots', runTimeParams.generatePlots);
    UnitTest.validationData('validationData2', validationData2);

end