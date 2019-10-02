function IBIOCDValidateFullAll(varargin)
% Perform a full data check of all validation functions
%
% Syntax:
%   IBIOCDValidateFullAll([varargin])
%
% Description:
%    Full data check (no figures, no publish) of all validation functions
%
%    This function contains examples of usage inline. To access these, type
%    'edit IBIOCFValidateFullAll' into the Command Window.
%
% Inputs:
%    None required.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    verbosity            - String. How chatty to be about the output.
%                           Default 'low'. The options are:
%           none: Don't say anything.
%           low: Minimal.
%           medium: As the name suggests.
%           high: More than medium.
%           max: As much as possible
%    generatePlots        - Boolean. Whether or not to generate the plots.
%                           Default false.
%    graphMismatchedData' - Boolean. Whether or not to makle a graph when
%                           the validation fails. Default true.
%    numericTolerance     - Numeric. The tolerance for numeric checks.
%                           Default 500 * eps.
%    asAssertion          - Boolean. Whether or not to throw an error if
%                           the validation fails. Default false.
%

% Examples:
%{
    IBIOColorDetectValidateFullAll('verbosity', 'high');
    IBIOColorDetectValidateFullAll('Numeric Tolerance', 1000 * eps);
    IBIOColorDetectValidateFullAll('generate plots', true);
%}

%% Close all figures so that we start with a clean slate
close all;

%% We will use preferences for the 'isetbioValidation' project
thisProject = 'IBIOColorDetect';
UnitTest.usePreferencesForProject(thisProject);

%% Parse input and set settable prefs
p = inputParser;
p.addParameter('verbosity', 'low', @ischar);
p.addParameter('generatePlots', false, @islogical);
p.addParameter('graphMismatchedData', false, @islogical);
p.addParameter('numericTolerance', 500 * eps, @isnumeric);
p.addParameter('asAssertion', false, @islogical);
p.parse(varargin{:});

UnitTest.setPref('verbosity', p.Results.verbosity);
UnitTest.setPref('generatePlots', p.Results.generatePlots);
UnitTest.setPref('graphMismatchedData', p.Results.graphMismatchedData);
UnitTest.setPref('numericTolerance', p.Results.numericTolerance);

%% Set other preferences for this function
% Run time error behavior
% valid options: 'rethrowExceptionAndAbort', 'catchExceptionAndContinue'
UnitTest.setPref('onRunTimeErrorBehavior', 'catchExceptionAndContinue');

% Plot generation
UnitTest.setPref('closeFigsOnInit', true);

%% Print current values of isetbioValidation prefs
UnitTest.listPrefs();

%% What to validate
listingScript = UnitTest.getPref('listingScript');
vScriptsList = eval(listingScript);

%% How to validate
% Run a FULL validation session (comparing actual data)
obj = UnitTest.runValidationSession(vScriptsList, 'FULLONLY');

if p.Results.asAssertion
    % assert no failed validations
    summary = [obj.summaryReport{:}];
    success = ~any([summary.fullFailed]);
    assert(success, 'One or more validations failed.');
end

end