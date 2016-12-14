function IBIOColorDetect()
%% Configure things for working on the IBIOColorDetect project.
%
% For use with the ToolboxToolbox.
%
% Copy this file into your ToolboxToolbox localToolboxHooks directory (by
% defalut, ~/localToolboxHooks) and rename it to "IBIOColorDetect.m".  Then
% edit your local copy in the section marked "Edit here", to match your
% local machine.
%
% Once you've copied and edited this file, you can set up the project with
%   tbUse('IBIOColorDetect')
%
% The thing that this does is add subfolders of the project to the Matlab
% path as well as define Matlab preferences that specify input and output
% directories.
%
% IBIOColorDetect()


%% Edit here.
%   By default, this will check environment variables and fall back on
%   something sensible.  Or, you can edit this section with local values.

% parent folder where IBIOColorDetect was cloned
projectBaseDir = getenv('IBIOCD_PROJECT_DIR');
if isempty(projectBaseDir)
    projectBaseDir = getpref('ToolboxToolbox', 'toolboxRoot');
end

% where to put generated output
outputBaseDir = getenv('IBIOCD_OUTPUT_DIR');
if isempty(outputBaseDir)
    outputBaseDir = fullfile(tbUserFolder(), 'IBIOColorDetect');
end


%% You should not need to edit below.


%% Say hello
fprintf('Running IBIOColorDetect local hook\n');
projectName = 'IBIOColorDetect';


%% Put IBIOColorDetect code onto the Matlab path
% Specify project URL and location
projectUrl = 'https://github.com/isetbio/IBIOColorDetect.git';

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
config = [withToolbox withTutorials withCompute withValidations];
tbDeployToolboxes('config', config, 'toolboxRoot', projectBaseDir, 'runLocalHooks', false);


%% Set up the UnitTestToolbox and RemoteDataToolbox.

remoteDataConfigFile = fullfile(projectBaseDir, ...
    'IBIOColorDetect', 'configuration', 'rdt-config-IBIOColorDetect.json');

p = struct(...
    'projectName',           projectName, ...                                                                                 % The project's name (also the preferences group name)
    'validationRootDir',     IBIOCDValidationDir(), ...                                                                       % Directory location where the 'scripts' subdirectory resides.
    'alternateFastDataDir',  '',  ...                                                                                         % Alternate FAST (hash) data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/fast
    'alternateFullDataDir',  '', ...                                                                                          % Alternate FULL data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/full
    'useRemoteDataToolbox',  true, ...                                                                                        % If true use Remote Data Toolbox to fetch full validation data on demand.
    'remoteDataToolboxConfig', remoteDataConfigFile, ...                                                                      % Struct, file path, or project name with Remote Data Toolbox configuration.
    'clonedWikiLocation',    '', ...                                                                                          % Local path to the directory where the wiki is cloned. Only relevant for publishing tutorials.
    'clonedGhPagesLocation', '', ...                                                                                          % Local path to the directory where the gh-pages repository is cloned. Only relevant for publishing tutorials.
    'githubRepoURL',         '', ...                                                                                          % Github URL for the project. This is only used for publishing tutorials.
    'generateGroundTruthDataIfNotFound',false,...                                                                             % Flag indicating whether to generate ground truth if one is not found
    'listingScript',         'IBIOCDValidateListAllValidationDirs', ...                                                       % Script that lists dirs to find validation scripts in
    'coreListingScript',     '', ...                                                                                          % Not used in this project
    'numericTolerance',      1e-11 ...                                                                                        % Numeric tolerance for comparisons with validation data.
    );

generatePreferenceGroup(p);
UnitTest.usePreferencesForProject(p.projectName);


%% Setup the output folder.
setpref(projectName, 'outputBaseDir', outputBaseDir);

end


%% Remove existing preferences and replace with given struct.
function generatePreferenceGroup(p)
if ispref(p.projectName)
    rmpref(p.projectName);
end

setpref(p.projectName, 'projectSpecificPreferences', p);
fprintf('Generated and saved preferences specific to the ''%s'' project.\n', p.projectName);
end
