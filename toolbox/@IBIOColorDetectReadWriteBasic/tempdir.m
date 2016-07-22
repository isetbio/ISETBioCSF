function tempDir = tempdir(obj,varagin)
% data = read(obj,name,varagin)
%
% Read data objects using metadata.

% Parse input.
p = inputParser;
p.parse(params,varargin{:});

% Find project directory using preferences
outputBaseDir = getpref('IBIOColorDetect','outputBaseDir');
tempDir = fullfile(outputBaseDir,'temp');
if (~exist(tempDir,'dir'))
    mkdir(tempDir);
end

end