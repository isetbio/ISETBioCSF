function IBIOColorDetectProjectConfig
% IBIOColorDetectConfig
%
% Configure things for working on the IBIOColorDetect project.

%% Setup Matlab environment using the ToolboxToolbox
tbDeployToolboxes('configPath','IBIOColorDetectToolboxConfig.json','resetPath',true);

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
