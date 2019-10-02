function varargout = v_IBIOCD2ColorDetectFindPerformance(varargin)
% Validate the ISETBioCSF color detect find performance tutorial.
%
% Syntax:
%   [varargout] = v_IBIOCDC2olorDetectFindPerformance([varargin])
%
% Description:
%    Works by running t_colorDetectFindPerformance with various arguments
%    and comparing results with those stored.
%
%    The 2 in the filename is to make sure that's gets run in the right
%    order on a run through all validation scripts.
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
        '***** v_IBIOCD2ColorDetectFindPerformance *****');

    %% Get run params
    rParams = responseParamsGenerate('fastComputeParams', true);

    %% Basic validation
    [validationData1, extraData1] = t_colorDetectFindPerformance(...
        'rParams', rParams, ...
        'employStandardHostComputerResources', true, ...
        'generatePlots', runTimeParams.generatePlots, ...
        'plotPsychometric', true, 'freezeNoise', true);
    UnitTest.validationData('validationData1', validationData1, ...
        'UsingTheFollowingVariableTolerancePairs', ...
         'validationData1.percentCorrect', 0.11, ...
         'validationData1.stdErr', 0.06);
    UnitTest.extraData('extraData1', extraData1);
end