function theDir = IBIOColorDetectOutputDir(conditionDir,subDir)
% theDir = IBIOColorDetectOutputDir(conditionDir,subDir)
%
% Return where output files should go.  This uses a base directory
% determined by the IBIOColorDetect preference 'outputBaseDir' if it
% exists, and otherwise creates a directory IBIOColorDetectOutput at the
% same level as IBIOColorDetect.
%
% See also IBIOColorDetectPreferencesTemplate.

if (nargin < 1 || isempty(conditionDir))
    conditionDir = '';
end
if (nargin < 2 || isempty(subDir))
    subDir = '';
end

if (ispref('IBIOColorDetect','outputBaseDir'))
    topDir = fullfile(getpref('IBIOColorDetect','outputBaseDir'),conditionDir);
else
    [p,~] = fileparts(which(mfilename()));
    topBaseDir = fullfile(p(1:strfind(p,'IBIOColorDetect')+numel('IBIOColorDetect')-1),conditionDir);
    topBaseDir = fullfile(topBaseDir,'..','IBIOColorDetectOutput');
    if (~exist(topBaseDir,'dir'))
        mkdir(topBaseDir);
    end
    topDir = fullfile(topBaseDir,conditionDir);
end

if (~exist(topDir,'dir'))
    mkdir(topDir);
end
theDir = fullfile(topDir,subDir);
if (~exist(theDir, 'dir'))
    mkdir(theDir);
end

end
