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
%   'FigureType'
%      Type of figure to save, as supported by FigureSave/saveas
%      'pdf' - PDF (default)
%   'MovieExtension'
%      Extension for movie filename
%      'm4v' - MPEG-4 (default)

%% Parse input.
p = inputParser;
p.addRequired('name',@ischar);
p.addRequired('data');
p.addParameter('Type','mat',@ischar);
p.addParameter('ArtifactParams','mat',@isstruct);
p.addParameter('FigureType','pdf',@ischar);
p.addParameter('MovieExtension','m4v',@ischar);
p.parse(name,data,varargin{:});

%% Get parent directory list and make sure the output tree is in place
theParentDir = '';
if (~isempty(obj.parentParamsList))
    for ii = 1:length(obj.parentParamsList)
        thisParentParams = obj.parentParamsList{ii};
        switch(obj.parentParams.type)
            case 'ResponseGeneration'
                thisParentDir = obj.paramsToResponseGenerationDirName(thisParentParams);
            case 'ColorModulation'
                thisParentDir = obj.paramsToColorModulationDirName(thisParentParams);
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
theCurrentDir = theParentDir;
for ii = 1:length(obj.currentParamsList)
    thisCurrentParams = obj.currentParamsList{ii};
    switch(thisCurrentParams.type)
        case 'ResponseGeneration'
            thisCurrentDir = obj.paramsToResponseGenerationDirName(thisCurrentParams);
        case 'ColorModulation'
            thisCurrentDir = obj.paramsToColorModulationDirName(thisCurrentParams);
        otherwise
            error('Unkown current parameters type');
    end
    theCurrentDir = fullfile(theCurrentDir,thisCurrentDir);
    if (~exist(theCurrentDir,'dir'))
        mkdir(theCurrentDir);
    end
end

%% Compose directory where data should go
switch (p.Results.Type)
    case 'mat'
        theDir = fullfile(theCurrentDir,'matfiles');
    case 'figure'
        theDir = fullfile(theCurrentDir,'figures');
    case {'movie','movieFile'}
        theDir = fullfile(theCurrentDir,'movies');
end
if (~exist(theDir,'dir'))
    mkdir(theDir);
end

%% Write the data
switch (p.Results.Type)
    case 'mat'
        save(fullfile(theDir,name),data,'-v7.3');
    case 'figure'
        if (exist('FigureSave','file'))
            FigureSave(fullfile(theDir,name),data,p.Results.FigureType);
        else
            saveas(h,fullfile(theDir,name),p.Results.FigureType);
        end
    case 'movieFile'
        unix(['mv ' data ' ' fullfile(theDir,name) '.' p.Results.MovieExtension]);
end

