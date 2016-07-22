function write(obj,name,data,varargin)
% write(obj,name,data,varargin)
%
% Write data objects in a good place, with metadata.
%
% Key/value pairs
%   'Type'
%     'mat' - Matlab .mat file
%     'figure' - A figure
%     'movie' - A movie
%     'movieFile' - String with full path to temp movie file
%   'ArtifactParams'
%      Structure with artifact specific information.

%% Parse input.
p = inputParser;
p.addRequired('name',@ischar);
p.addRequired('data');
p.addParameter('Type','mat',@ischar);
p.addParameter('ArtifactParams','mat',@isstruct);
p.parse(name,data,varargin{:});

%% Get parent directory string
if (~isempty(obj.parentParams))
    switch(obj.parentParams.type)
        case 'ResponseGeneration'
            parentDir = responseGenerationDirName(obj.parentParams);
        otherwise
            error('Unkown parent parameters type');
    end
else
    parentDir = '';
end

%% Get current dir
switch(obj.currentParams.type)
    case 'ResponseGeneration'
        currentDir = responseGenerationDirName(obj.parentParams);
    otherwise
        error('Unkown parent parameters type');
end

%% Compose where data should go
theDir = fullfile(getpref('IBIOColorDetect','outputBaseDir'),parentDir,currentDir);

switch (p.Results.Type)
    case 'movieFile'
end

videoDir = colorGaborDetectOutputDir(conditionDir,'videos');

end