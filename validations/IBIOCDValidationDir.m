function validationDir = IBIOCDValidationDir
% Return the full ISETBioCSF validation directory path
%
% Syntax:
%   validationDir = IBIOCDValidationDir()
%
% Description:
%    Return the full path to the IBIOColorDetect project validtion dir,
%    which should be where this function resides in the toolbox.
%
% Inputs:
%    None.
%
% Outputs:
%    validationDir - String. The full path to the validation directory.
%
% Optional key/value pairs:
%    None.
%

fullValidationDir = mfilename('fullpath');
validationDir = fileparts(fullValidationDir);

end
