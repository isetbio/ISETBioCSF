function IBIOCDRunTutorialsAll
% Run all of the linked tutorials
%
% Syntax:
%   IBIOCDRunTutorialsAll
%
% Description:
%    Run all of the ISETBioCSF tutorials that we think should work, and
%    print out a report at the end as to whether they threw errors, or not.
%    Scripts inside of RootPath/tutorials are run, except that scripts
%    within the directory 'underDevelopment' or 'support' are skipped.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    07/26/17  dhb  Wrote this, because we care.
%    10/01/19  jnm  Documentation pass

%% User/project specific preferences
% tutorialsSourceDir - A local directory where tutorial scripts are located
p = struct('rootDirectory', fileparts(which(mfilename())), ...
    'tutorialsSourceDir', ...
    fullfile(fileparts(which(mfilename())), '..', 'tutorials'));

%% List of scripts to be skipped from automatic publishing.
% Anything with this in its path name is skipped.
scriptsToSkip = {'underDevelopment', 'support'};

%% Use UnitTestToolbox method to do this.
UnitTest.runProjectTutorials(p, scriptsToSkip, 'All');

end
