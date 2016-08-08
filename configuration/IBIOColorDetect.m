function IBIOColorDetect
% IBIOColorDetect
%
% Configure things for working on the IBIOColorDetect project.

%% Say hello
fprintf('Running IBIOColorDetect local hook\n');

%% Put project toolbox onto path
tbDeployToolboxes('config',tbToolboxRecord( ...
    'name', 'IBIOColorDetectProjectToolbox', ...
    'type', 'local', ...
    'url', fullfile(pwd,'toolbox')) ...
    );

%% Set preferences for project output
%
% Root dir
%baseDir = '/Users/dhb/DropboxLab/IBIO_analysis';
baseDir = '/Volumes/Users1/DropboxLab/IBIO_analysis';
if (~exist(baseDir))
    mkdir(baseDir);
end

% This project's dir under the base dir
projectDir = 'IBIOColorDetect';
theDir = fullfile(baseDir,projectDir);

% Set the preference
setpref('IBIOColorDetect','outputBaseDir',theDir);

%% Run isetbio prefs
run('/Users/dhb/Desktop/ProjectPrefs/isetbio/ieSetPreferencesDB');
