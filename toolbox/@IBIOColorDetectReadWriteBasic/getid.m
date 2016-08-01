function fileid = getid(obj,name,parentParamsList,currentParamsList,theProgram,varargin)
% fileid = getid(obj,name,parentParamsList,currentParamsList,theProgram,varargin)
%
% Get a unique id required for reading and writing a piece of data.
% Here that is the full path to the file, but in a database world it could
% be something else.
%
% Key/value pairs
%   'Type'
%     'mat' - Matlab .mat file
%     'figure' - A figure
%     'movie' - A movie
%     'movieFile' - String with full path to temp movie file
%   'ArtifactParams'
%      Structure with artifact specific information.
%   'FigureType'
%      Type of figure to save, represented as file extension.
%      'pdf' - PDF (default)
%   'MovieType'
%      Type for movie to save, represented as file extension.
%      'm4v' - MPEG-4 (default)

%% Parse input.
p = inputParser;
p.addRequired('name',@ischar);
p.addRequired('parentParamsList',@iscell);
p.addRequired('currentParamsList',@iscell);
p.addRequired('theProgram',@ischar);
p.addParameter('Type','mat',@ischar);
p.addParameter('ArtifactParams',[],@isstruct);
p.addParameter('FigureType','pdf',@ischar);
p.addParameter('MovieType','m4v',@ischar);
p.parse(name,parentParamsList,currentParamsList,theProgram,varargin{:});

%% Get parent directory list and make sure the output tree is in place
theParentDir = fullfile(getpref('IBIOColorDetect','outputBaseDir'));
if (~exist(theParentDir,'dir'))
    mkdir(theParentDir);
end
if (~isempty(p.Results.parentParamsList))
    for ii = 1:length(p.Results.parentParamsList)
        thisParentParams = p.Results.parentParamsList{ii};
        switch(thisParentParams.type)
            case 'ResponseGeneration'
                thisParentDir = obj.paramsToResponseGenerationDirName(thisParentParams);
            case 'ColorModulation'
                thisParentDir = obj.paramsToColorModulationDirName(thisParentParams);
            case 'LMPlaneInstance'
                thisParentDir = obj.paramsToLMPlaneInstanceDirName(thisParentParams);
            case 'threshold'
                thisParentDir = obj.paramsToThresholdDirName(thisParentParams);
            case 'psychoEllipsoid'
                thisParentDir = obj.paramsToPsychoEllipsoidDirName(thisParentParams);
            otherwise
                error('Unkown parent parameters type');
        end
        theParentDir = fullfile(theParentDir,thisParentDir);
        if (~exist(theParentDir,'dir'))
            mkdir(theParentDir);
        end
    end    
end

%% Get current dir name and make it if necessary
theCurrentDir = theParentDir;
for ii = 1:length(p.Results.currentParamsList)
    thisCurrentParams = p.Results.currentParamsList{ii};
    switch(thisCurrentParams.type)
        case 'ResponseGeneration'
            thisCurrentDir = obj.paramsToResponseGenerationDirName(thisCurrentParams);
        case 'ColorModulation'
            thisCurrentDir = obj.paramsToColorModulationDirName(thisCurrentParams);
        case 'LMPlaneInstance'
            thisCurrentDir = obj.paramsToLMPlaneInstanceDirName(thisCurrentParams);
        case 'threshold'
            thisCurrentDir = obj.paramsToThresholdDirName(thisCurrentParams);
        case 'psychoEllipsoid'
            thisCurrentDir = obj.paramsToPsychoEllipsoidDirName(thisCurrentParams);
        otherwise
            error('Unkown current parameters type');
    end
    theCurrentDir = fullfile(theCurrentDir,thisCurrentDir);
    if (~exist(theCurrentDir,'dir'))
        mkdir(theCurrentDir);
    end
end

%% Filetype specific subdirectory
switch (p.Results.Type)
    case 'mat'
        theTypeDir = fullfile(theCurrentDir,'matfiles');
    case 'figure'
        theTypeDir = fullfile(theCurrentDir,'figures');
    case {'movie','movieFile'}
        theTypeDir = fullfile(theCurrentDir,'movies');
end
if (~exist(theTypeDir,'dir'))
    mkdir(theTypeDir);
end

%% Add in the name of reading or writing program
theDir = fullfile(theTypeDir,theProgram);
if (~exist(theDir,'dir'))
    mkdir(theDir);
end

%% Add in the filename
switch (p.Results.Type)
    case 'mat'
        fileid = fullfile(theDir,[name '.mat']);
    case 'figure'
        fileid = fullfile(theDir,[name '.' p.Results.FigureType]);
    case {'movie','movieFile'}
        fileid = fullfile(theDir,[name '.' p.Results.MovieType]);
end

end