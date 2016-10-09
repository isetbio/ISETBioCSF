function validationDir = tfeValidationDir
%
% Return the full path to the tfe project validtion dir, which should be
% where this function resides in the toolbox.

fullValidationDir = mfilename('fullpath');
validationDir = fileparts(fullValidationDir);

end
