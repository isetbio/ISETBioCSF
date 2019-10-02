function vScriptsList = IBIOCDValidateListAllValidationDirs
% List specified validation directories
%
% Syntax:
%   vScriptsList = IBIOCDValidateListAllValidationDirs
%
% Description:
%    This encapsulates a full list of our validation directories, so we
%    only need to update it in one place.
%
%    Note: This doesn't list the example scripts, and doesn't override any
%    default prefs.
%
%    List of script directories to validate. Each entry contains a cell
%    array with with a validation script directory and an optional struct
%    with prefs that override the corresponding isetbioValidation prefs.
%    At the moment only the 'generatePlots' pref can be overriden.
%
% Inputs:
%    None.
%
% Outputs:
%    vScriptsList - Cell. A cell list of the validation script directories.
%                   Currently only includes basic & compute.
%
% Optional key/value pairs:
%    None.
%

% Get rootDir
rootDir = UnitTest.getPref('validationRootDir');

vScriptsList = {{fullfile(rootDir, 'scripts', 'basic')} ...
    {fullfile(rootDir, 'scripts', 'compute')}};

end
