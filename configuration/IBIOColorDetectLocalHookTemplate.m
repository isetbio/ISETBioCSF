function IBIOColorDetectLocalHook
% Configure things for working on the IBIOColorDetect project.
%
% Syntax:
%   IBIOColorDetectLocalHook
%
% Description:
%    For use with the ToolboxToolbox.
%
%    If you 'git clone' IBIOColorDetect into your ToolboxToolbox
%    "projectRoot" folder, then run in MATLAB
%
%        tbUseProject('IBIOColorDetect')
%
%    ToolboxToolbox will set up IBIOColorDetect and its dependencies on
%    your machine.
%
%    As part of the setup process, ToolboxToolbox will copy this file to your
%    ToolboxToolbox localToolboxHooks directory (minus the "Template" suffix).
%    The defalt location for this would be
%
%        ~/localToolboxHooks/IBIOColorDetectLocalHook.m
%
%    Each time you run tbUseProject('IBIOColorDetect'), ToolboxToolbox will
%    execute your local copy of this file to do setup for IBIOColorDetect.
%
%    You should edit your local copy with values that are correct for your
%    local machine, for example the output directory location.
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

%% Say hello.
fprintf('IBIOColorDetect local hook.\n');
projectName = 'IBIOColorDetect';

%% UnitTestToolbox and RemoteDataToolbox setup.
% If you customize your rdt-config json file, you will want to place it
% somewhere outside the repository and change the path to point to your
% copy.  The usual reason for doing this is so that you can add your
% username and password and have write permission.
projectBaseDir = tbLocateProject(projectName);
rdtConfig = fullfile(projectBaseDir, 'configuration', ...
    ['rdt-config-' projectName '.json']);
% Define the project parameters inside the structure p.
%    projectName: The project's name (also the preferences group name)
%    validationRootDir: The directory location where the 'scripts'
%                       subdirectory resides.
%    alternateFastDataDir: An alternate FAST (hash) data directory
%                          location. Specify '' to use the default
%                          location, i.e., $validationRootDir/data/fast
%    alternateFullDataDir: An alternate FULL data directory location.
%                          Specify '' to use the default location, i.e.,
%                          $validationRootDir/data/full
%    useRemoteDataToolbox: A boolean, which if true uses Remote Data
%                          Toolbox to fetch full validation data on demand.
%    remoteDataToolboxConfig: A struct, file path, or project name with
%                             Remote Data Toolbox configuration.
%    clonedWikiLocation: A local path to the directory where the wiki is
%                        cloned. Only relevant for publishing tutorials.
%    clonedGhPagesLocation: A local path to the directory where the
%                           gh-pages repository is cloned. Only relevant
%                           for publishing tutorials.
%    githubRepoURL: The Github URL for the project. This is only used for
%                   publishing tutorials.
%    generateGroundTruthDataIfNotFound: The boolean flag indicating whether
%                                       to generate ground truth if one is
%                                       not found.
%    listingScript: The script that lists the directory in which to find
%                   the validation scripts.
%    coreListingScript: This is not used in this project
%    numericTolerance: The numeric tolerance for comparisons with
%                      validation data.
p = struct(...
    'projectName', projectName, ...
    'validationRootDir', IBIOCDValidationDir, ...
    'alternateFastDataDir', '', ...
    'alternateFullDataDir', '', ...
    'useRemoteDataToolbox', true, ...
    'remoteDataToolboxConfig', rdtConfig, ...
    'clonedWikiLocation', '', ...
    'clonedGhPagesLocation', '', ...
    'githubRepoURL', '', ...
    'generateGroundTruthDataIfNotFound', true, ...
    'listingScript', 'IBIOCDValidateListAllValidationDirs', ...
    'coreListingScript', '', ...
    'numericTolerance', 1e-11);

generatePreferenceGroup(p);
UnitTest.usePreferencesForProject(p.projectName);

%% Output directory.
% This is where the project writes its output. By default, we'll stick it
% in a subfolder of a folder called output, in the tbUserFolder. But you
% may want it somewhere else.
outputBaseDir = fullfile(tbUserFolder(), 'output', projectName);
if (7 ~= exist(outputBaseDir, 'dir')), mkdir(outputBaseDir); end
setpref(projectName, 'outputBaseDir', outputBaseDir);

%% Backwards compatibility prefs
% setpref('IBIOColorDetectBackCompat','oiWithCustomOptics',true);

end

function generatePreferenceGroup(p)
% Generate preferences that work with UnitTest toolbox.
%
% Syntax:
%   generatePreferenceGroup(p)
%
% Description:
%    Generate a set of preferences that will work with UnitTest toolbox.
%
% Inputs:
%    p - Struct. A structure of preferences.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% Remove any existing preferences for this project
if ispref(p.projectName), rmpref(p.projectName); end

% Renerate and save the project-specific preferences
setpref(p.projectName, 'projectSpecificPreferences', p);
fprintf(['Generated and saved preferences specific to the ''%s''' ...
    ' project.\n'], p.projectName);
end
