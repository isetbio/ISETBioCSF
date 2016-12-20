function IBIOCDValidateFullOne(varargin)
% IBIOCDValidateFullOne(varargin)
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

%% Close all figures so that we start with a clean slate
close all; 

%% We will use preferences for the 'isetbioValidation' project
thisProject = 'IBIOColorDetect';
UnitTest.usePreferencesForProject(thisProject);

%% Parse input and set settable prefs
p = inputParser;
p.addParameter('verbosity','high',@ischar);
p.addParameter('generatePlots',false,@islogical);
p.addParameter('graphMismatchedData',true,@islogical);
p.addParameter('numericTolerance',500*eps,@isnumeric);
p.parse(varargin{:});
UnitTest.setPref('verbosity',p.Results.verbosity);
UnitTest.setPref('generatePlots',p.Results.generatePlots);
UnitTest.setPref('graphMismatchedData',p.Results.graphMismatchedData);
UnitTest.setPref('numericTolerance',p.Results.numericTolerance);

% Run time error behavior
% valid options are: 'rethrowExceptionAndAbort', 'catchExceptionAndContinue'
UnitTest.setPref('onRunTimeErrorBehavior', 'rethrowExceptionAndAbort');

% Plot generation
UnitTest.setPref('closeFigsOnInit', true);

%% Print all existing validation scripts and ask the user to select one for validation
singleScriptToValidate = UnitTest.selectScriptFromExistingOnes();

%% Validate
UnitTest.runValidationSession({{singleScriptToValidate, []}}, 'FULLONLY');

end