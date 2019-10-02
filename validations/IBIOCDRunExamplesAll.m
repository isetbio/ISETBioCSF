%% IBIOCDRunExamplesAll
% Run all the examples in the ISETBioCSF tree
%
% Syntax:
%   IBIOCDRunExamplesAll
%
% Description:
%     Run all the examples in the ISETBioCSF tree, excepthose that contain
%     a line of the form "% ETTBSkip"
%
% See Also:
%   ieValidateFullAll, IBIOCDRunTutorialsAll

% History:
%    01/17/18  dhb  Wrote it.
%    10/01/19  jnm  Documentation pass

ExecuteExamplesInDirectory(fullfile(...
    fileparts(which(mfilename())), '..'), 'verbose', false);
