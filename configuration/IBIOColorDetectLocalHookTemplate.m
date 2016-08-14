function IBIOColorDetect
% IBIOColorDetect
%
% Configure things for working on the IBIOColorDetect project.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by defalut,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'AOMontaging'}) to set up for
% this project.  You then edit your local copy to match your 
%
% The thing that this does is add subfolders of the project to the path as
% well as define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

%% Say hello
fprintf('Running IBIOColorDetect local hook\n');

%% Put project toolbox onto path
%
% Specify project name and location
projectName = 'IBIOColorDetect';
projectBaseDir = '/Users/Shared/Matlab/Analysis';

% Put the subdirs of the project we need on the path
tbDeployToolboxes('config',tbToolboxRecord( ...
    'name', 'IBIOColorDetectProjectToolbox', ...
    'type', 'local', ...
    'url', fullfile(projectBaseDir,projectName,'toolbox')) ...
    );
tbDeployToolboxes('config',tbToolboxRecord( ...
    'name', 'IBIOColorDetectProjectTutorials', ...
    'type', 'local', ...
    'url', fullfile(projectBaseDir,projectName,'tutorials')) ...
    );

%% Set preferences for project output
%
%outputBaseDir = '/Users/dhb/DropboxLab/IBIO_analysis';
outputBaseDir = '/Volumes/Users1/DropboxLab/IBIO_analysis';

% Make base directory if it doesn't exist
if (~exist(outputBaseDir))
    mkdir(outputBaseDir);
end

% This project's dir under the base dir
theDir = fullfile(outputBaseDir,projectName);

% Set the preference
setpref('IBIOColorDetect','outputBaseDir',theDir);
