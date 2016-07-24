function tempdir = tempdir(obj,varargin)
% tempdir = tempdir(obj,varargin)
%
% Read data objects using metadata.

% Parse input.
p = inputParser;
p.parse(varargin{:});

% Find project directory using preferences
outputBaseDir = getpref('IBIOColorDetect','outputBaseDir');
tempdir = fullfile(outputBaseDir,'temp');
if (~exist(tempdir,'dir'))
    mkdir(tempdir);
end

end