function IBIOCDValidateFullAll(varargin)
% IBIOCDValidateFullAll(varargin)
%
% Full data check (no figures, no publish) of all validation functions
%
% Optional key/value pairs
%   'verbosity' - string (default 'low').  How chatty to be about output.
%      'none' - Don't say anything.
%      'low' - Minimal.
%      'medium' - As the name suggests.
%      'high' - More than medium.
%      'max' - As much as possible
%   'generatePlots' - true/false (default false).  Generate plots?
%   'graphMismatchedData' - true/false (default true).  Make a graph when
%       validation fails?
%   'numericTolerance' - value (default 500*eps).  Tolerance to use for numeric checks.
%   'asAssertion' - true/false (default false).  Throw error if validation fails?

% Examples:
%   IBIOColorDetectValidateFullAll('verbosity','high');
%   IBIOColorDetectValidateFullAll('Numeric Tolerance',1000*eps);
%   IBIOColorDetectValidateFullAll('generate plots',true);

%% Close all figures so that we start with a clean slate
close all;

assert(false, 'Testing validation failure.');

%% We will use preferences for the 'isetbioValidation' project
thisProject = 'IBIOColorDetect';
UnitTest.usePreferencesForProject(thisProject);

%% Parse input and set settable prefs
p = inputParser;
p.addParameter('verbosity','low',@ischar);
p.addParameter('generatePlots',false,@islogical);
p.addParameter('graphMismatchedData',false,@islogical);
p.addParameter('numericTolerance',500*eps,@isnumeric);
p.addParameter('asAssertion',false,@islogical);
p.parse(varargin{:});
UnitTest.setPref('verbosity',p.Results.verbosity);
UnitTest.setPref('generatePlots',p.Results.generatePlots);
UnitTest.setPref('graphMismatchedData',p.Results.graphMismatchedData);
UnitTest.setPref('numericTolerance',p.Results.numericTolerance);

%% Set other preferences for this function

% Run time error behavior
% valid options are: 'rethrowExceptionAndAbort', 'catchExceptionAndContinue'
UnitTest.setPref('onRunTimeErrorBehavior', 'catchExceptionAndContinue');

% Plot generation
UnitTest.setPref('closeFigsOnInit', true);

%% Print current values of isetbioValidation prefs
UnitTest.listPrefs();

%% What to validate
listingScript = UnitTest.getPref('listingScript');
vScriptsList = eval(listingScript);

%% How to validate
%
% Run a FULL validation session (comparing actual data)
obj = UnitTest.runValidationSession(vScriptsList, 'FULLONLY');

if p.Results.asAssertion
    % assert no failed validations
    summary = [obj.summaryReport{:}];
    success = ~any([summary.fullFailed]);
    assert(success, 'One or more validations failed.');
end

end