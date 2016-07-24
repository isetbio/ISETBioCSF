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
%   'MovieExtension'
%      Extension for movie filename
%      'm4v' - MPEG-4 (default)

%% Parse input.
p = inputParser;
p.addRequired('name',@ischar);
p.addRequired('data');
p.addParameter('Type','mat',@ischar);
p.addParameter('ArtifactParams','mat',@isstruct);
p.addParameter('MovieExtension','m4v',@isstruct);
p.parse(name,data,varargin{:});

%% Get parent directory list and make sure the output tree is in place
theParentDir = '';
if (~isempty(obj.parentParamsList))
    for ii = 1:length(obj.parentParamsList)
        thisParentParams = obj.parentParamsList{ii};
        switch(obj.parentParams.type)
            case 'ResponseGeneration'
                thisParentDir = obj.paramsToResponseGenerationDirName(thisParentParams);
            otherwise
                error('Unkown parent parameters type');
        end
        theParentDir = fullfile(theParentDir,thisParentDir);
        if (~exist(theParentDir,'dir'))
            mkdir(theParentDir);
        end
    end
    theParentDir = fullfile(getpref('IBIOColorDetect','outputBaseDir'),theParentDir);
else
    theParentDir = getpref('IBIOColorDetect','outputBaseDir');
    if (~exist(theParentDir,'dir'))
        mkdir(theParentDir);
    end
end

%% Get current dir name and make it if necessary
switch(obj.currentParams.type)
    case 'ResponseGeneration'
        currentDir = obj.paramsToResponseGenerationDirName(obj.currentParams);
    otherwise
        error('Unkown parent parameters type');
end
theCurrentDir = fullfile(theParentDir,currentDir);
if (~exist(theCurrentDir,'dir'))
    mkdir(theCurrentDir);
end

%% Compose directory where data should go
switch (p.Results.Type)
    case 'mat'
        theDir = fullfile(theCurrentDir,'matfiles');
    case 'movieFile'
        theDir = fullfile(theCurrentDir,'movies');
end
if (~exist(theDir,'dir'))
    mkdir(theDir);
end

%% Write the data
switch (p.Results.Type)
    case 'mat'
        save(fullfile(theDir,name),data,'-v7.3');
    case 'movieFile'
        unix(['mv ' data ' ' fullfile(theDir,name) '.' p.Results.MovieExtension]);
end

