function IBIOColorDetectConfig
% IBIOColorDetectConfig
%
% Configure things for working on the IBIOColorDetect project.

%% Add project toolbox to Matlab path
AddToMatlabPathDynamically(fullfile(fileparts(which(mfilename)),'toolbox')); 

%% Set preferences

% Root dir
baseDir = '/Volumes/Users1/DropboxLab/IBIO_analysis';
if (~exist(baseDir))
    mkdir(baseDir);
end

% This project's dir
projectDir = 'IBIOColorDetect';
theDir = fullfile(baseDir,projectDir);

% Set preferences
setpref('IBIOColorDetect','outputBaseDir',theDir);