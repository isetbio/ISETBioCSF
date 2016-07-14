function IBIOColorDetectPreferencesTemplate
% 
%
% Set a preference for where this project should dump output

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

end

