function IBIOColorDetectLocalHookTemplate
% IBIOColorDetect
%
% Configure things for working on the IBIOColorDetect project.
%
% For use with the ToolboxToolbox.
%
% If you 'git clone' IBIOColorDetect into your ToolboxToolbox "projectRoot"
% folder, then run in MATLAB
%   tbUseProject('IBIOColorDetect')
% ToolboxToolbox will set up IBIOColorDetect and its dependencies on
% your machine.
%
% As part of the setup process, ToolboxToolbox will copy this file to your
% ToolboxToolbox localToolboxHooks directory (minus the "Template" suffix).
% The defalt location for this would be
%   ~/localToolboxHooks/IBIOColorDetectLocalHook.m
%
% Each time you run tbUseProject('IBIOColorDetect'), ToolboxToolbox will
% execute your local copy of this file to do setup for IBIOColorDetect.
%
% You should edit your local copy with values that are correct for your
% local machine, for example the output directory location.
%


%% Say hello.
fprintf('IBIOColorDetect local hook.\n');
projectName = 'IBIOColorDetect';


%% UnitTestToolbox and RemoteDataToolbox setup.
projectBaseDir = tbLocateProject(projectName);
rdtConfig = fullfile(projectBaseDir, 'configuration', ['rdt-config-' projectName '.json']);
p = struct(...
    'projectName',           projectName, ...                                                                                 % The project's name (also the preferences group name)
    'validationRootDir',     IBIOCDValidationDir, ...                                                                         % Directory location where the 'scripts' subdirectory resides.
    'alternateFastDataDir',  '',  ...                                                                                         % Alternate FAST (hash) data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/fast
    'alternateFullDataDir',  '', ...                                                                                          % Alternate FULL data directory location. Specify '' to use the default location, i.e., $validationRootDir/data/full
    'useRemoteDataToolbox',  true, ...                                                                                        % If true use Remote Data Toolbox to fetch full validation data on demand.
    'remoteDataToolboxConfig', rdtConfig, ...                                                                                 % Struct, file path, or project name with Remote Data Toolbox configuration.
    'clonedWikiLocation',    '', ...                                                                                          % Local path to the directory where the wiki is cloned. Only relevant for publishing tutorials.
    'clonedGhPagesLocation', '', ...                                                                                          % Local path to the directory where the gh-pages repository is cloned. Only relevant for publishing tutorials.
    'githubRepoURL',         '', ...                                                                                          % Github URL for the project. This is only used for publishing tutorials.
    'generateGroundTruthDataIfNotFound',true,...                                                                              % Flag indicating whether to generate ground truth if one is not found
    'listingScript',         'IBIOCDValidateListAllValidationDirs', ...                                                       % Script that lists dirs to find validation scripts in
    'coreListingScript',     '', ...                                                                                          % Not used in this project
    'numericTolerance',      1e-11 ...                                                                                        % Numeric tolerance for comparisons with validation data.
    );

generatePreferenceGroup(p);
UnitTest.usePreferencesForProject(p.projectName);


%% Output directory.
outputBaseDir = fullfile(tbUserFolder(), 'output', projectName);
if (7 ~= exist(outputBaseDir, 'dir'))
    mkdir(outputBaseDir);
end

setpref(projectName, 'outputBaseDir', outputBaseDir);

end


%% Generate preferences that work with UnitTest toolbox.
function generatePreferenceGroup(p)
% Remove any existing preferences for this project
if ispref(p.projectName)
    rmpref(p.projectName);
end

% Renerate and save the project-specific preferences
setpref(p.projectName, 'projectSpecificPreferences', p);
fprintf('Generated and saved preferences specific to the ''%s'' project.\n', p.projectName);
end
