function IBIOColorDetect
% IBIOColorDetect
%
% Configure things for working on the IBIOColorDetect project.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by defalut,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'IBIOColorDetect'}) to set up for
% this project.  You then edit your local copy to match your local machine.
%
% The thing that this does is add subfolders of the project to the path as
% well as define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

%% Say hello
fprintf('Running IBIOColorDetect local hook\n');
projectName = 'IBIOColorDetect';

%% Put project toolbox onto path
%
% Specify project URL and location
projectUrl = 'https://github.com/isetbio/IBIOColorDetect.git';
projectBaseDir = '/Users/Shared/Matlab/Analysis';

% Declare the project git repo and two subfolders that we want on the path
% Only update repository the first time, to save time.
withToolbox = tbToolboxRecord( ...
    'name', 'IBIOColorDetect', ...
    'type', 'git', ...
    'url', projectUrl, ...
    'subfolder', 'toolbox');
withTutorials = tbToolboxRecord( ...
    'name', 'IBIOColorDetect', ...
    'type', 'git', ...
    'url', projectUrl, ...
    'subfolder', 'tutorials', ...
    'update', 'never');
withCompute = tbToolboxRecord( ...
    'name', 'IBIOColorDetect', ...
    'type', 'git', ...
    'url', projectUrl, ...
    'subfolder', 'compute', ...
    'update', 'never');
withValidations = tbToolboxRecord( ...
    'name', 'IBIOColorDetect', ...
    'type', 'git', ...
    'url', projectUrl, ...
    'subfolder', 'validations', ...
    'update', 'never');

% Obtain or update the git repo and add subfolders to the Matlab path
config = [withToolbox withTutorials with Compute withValidations];
tbDeployToolboxes('config', config, 'toolboxRoot', projectBaseDir, 'runLocalHooks', false);

%% Specify project-specific preferences
%
% This currently include UnitTestToolbox/RemoteDataToolbox setup
p = struct(...
    'projectName',           projectName, ...                                                                                 % The project's name (also the preferences group name)
    'validationRootDir',     IBIOColorDetectValidationDir, ...                                                                % Directory location where the 'scripts' subdirectory resides.
    'alternateFastDataDir',  '',  ...                                                                                         % Alternate FAST (hash) data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/fast
    'alternateFullDataDir',  '', ...                                                                                          % Alternate FULL data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/full
    'useRemoteDataToolbox',  true, ...                                                                                        % If true use Remote Data Toolbox to fetch full validation data on demand.
    'remoteDataToolboxConfig', projectName, ...                                                                               % Struct, file path, or project name with Remote Data Toolbox configuration.
    'clonedWikiLocation',    '', ...                                                                                          % Local path to the directory where the wiki is cloned. Only relevant for publishing tutorials.
    'clonedGhPagesLocation', '', ...                                                                                          % Local path to the directory where the gh-pages repository is cloned. Only relevant for publishing tutorials.
    'githubRepoURL',         '', ...                                                                                          % Github URL for the project. This is only used for publishing tutorials.
    'generateGroundTruthDataIfNotFound',false,...                                                                             % Flag indicating whether to generate ground truth if one is not found
    'listingScript',         'IBIOColorDetectValidateListAllValidationDirs', ...                                              % Script that lists dirs to find validation scripts in
    'coreListingScript',     '', ...                                                                                          % Not used in this project
    'numericTolerance',      1e-11 ...                                                                                        % Numeric tolerance for comparisons with validation data.
    );

generatePreferenceGroup(p);
UnitTest.usePreferencesForProject(p.projectName);

%% Set up additional preferences for project output
%
%outputBaseDir = '/Users/dhb/DropboxLab/IBIO_analysis';
outputBaseDir = '/Volumes/Users1/DropboxLab/IBIO_analysis';
if (~exist(outputBaseDir))
    mkdir(outputBaseDir);
end

% This project's dir under the base dir
theDir = fullfile(outputBaseDir,projectName);

% Add into the preferences
setpref(projectName,'outputBaseDir',theDir);

end

%% Little generatePreferenceGroup utility
function generatePreferenceGroup(p)
% Remove any existing preferences for this project
if ispref(p.projectName)
    rmpref(p.projectName);
end

% Renerate and save the project-specific preferences
setpref(p.projectName, 'projectSpecificPreferences', p);
fprintf('Generated and saved preferences specific to the ''%s'' project.\n', p.projectName);
end
