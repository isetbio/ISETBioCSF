function IBIOCDValidateFullOne(varargin)
% Validate a single script
%
% Syntax:
%   IBIOCDValidateFullOne([varargin])
%
% Description:
%    Validate the selected script.
%
% Inputs:
%    None required.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    verbosity           - String. The level of verbosity for the output.
%                          Default 'high'. Options are:
%           none: Don't say anything.
%           low: Minimal.
%           medium: As the name suggests.
%           high: More than medium.
%           max: As much as possible
%    generatePlots       - Boolean. Whether to generate the requisite
%                          plots. Default false.
%    graphMismatchedData - Boolean. Whether to create a graph when the
%                          validation fails. Default true.
%    numericTolerance    - Numeric. The tolerance to use for numeric
%                          checks. Default 500 * eps.
%

%% Close all figures so that we start with a clean slate
close all;

%% We will use preferences for the 'IBIOColorDetect' project
thisProject = 'IBIOColorDetect';
UnitTest.usePreferencesForProject(thisProject);

%% Parse input and set settable prefs
p = inputParser;
p.addParameter('verbosity', 'high', @ischar);
p.addParameter('generatePlots', false, @islogical);
p.addParameter('graphMismatchedData', true, @islogical);
p.addParameter('numericTolerance', 500 * eps, @isnumeric);
p.parse(varargin{:});
UnitTest.setPref('verbosity', p.Results.verbosity);
UnitTest.setPref('generatePlots', p.Results.generatePlots);
UnitTest.setPref('graphMismatchedData', p.Results.graphMismatchedData);
UnitTest.setPref('numericTolerance', p.Results.numericTolerance);

% Run time error behavior
% valid options: 'rethrowExceptionAndAbort', 'catchExceptionAndContinue'
UnitTest.setPref('onRunTimeErrorBehavior', 'rethrowExceptionAndAbort');

% Plot generation
UnitTest.setPref('closeFigsOnInit', true);

%% Instigate user selection
% Print all of the existing validation scripts and ask the user to select
% one for validation.
singleScriptToValidate = UnitTest.selectScriptFromExistingOnes();

%% Validate
UnitTest.runValidationSession({{singleScriptToValidate, []}}, 'FULLONLY');

end
