function [fileid,filedir,filename] = getid(obj,name,paramsList,theProgram,varargin)
% [fileid,filedir,filename] = getid(obj,name,paramsList,theProgram,varargin)
%
% Get a unique id required for reading and writing a piece of data.
% Here that is the full path to the file, but in a database world it could
% be something else.  This also returns the directory and filename that
% make up the fileid.
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
%   'MakeDirectories' - true/false (default true).  Make the directories if
%      they don't exist.  This is very useful for writes, but not so good
%      for reads or deletes as it creates orphan directories if you don't
%      pass in the right arguments.

%% Parse input.
p = inputParser;
p.addRequired('name',@ischar);
p.addRequired('paramsList',@iscell);
p.addRequired('theProgram',@ischar);
p.addParameter('Type','mat',@ischar);
p.addParameter('ArtifactParams',[],@isstruct);
p.addParameter('FigureType','pdf',@ischar);
p.addParameter('FigureHandle', []);
p.addParameter('MovieType','m4v',@ischar);
p.addParameter('MakeDirectories',true,@islogical);
p.addParameter('ShowFilePath', false,@islogical);
p.parse(name,paramsList,theProgram,varargin{:});

%% Get parent directory list and make sure the output tree is in place
theParentDir = fullfile(getpref('IBIOColorDetect','outputBaseDir'));
if (~exist(theParentDir,'dir') & p.Results.MakeDirectories)
    mkdir(theParentDir);
end
if (~isempty(p.Results.paramsList))
    for ii = 1:length(p.Results.paramsList)
        thisParams = p.Results.paramsList{ii};
        switch(thisParams.type)
            case 'TopLevelDir'
                thisParentDir = obj.paramsToTopLevelDirName(thisParams);
                
            case 'ResponseGeneration'
                thisParentDir = obj.paramsToResponseGenerationDirName(thisParams);
                
            case {'ColorModulation'}
                thisParentDir = obj.paramsToColorModulationDirName(thisParams);
                
            case {'Background'}
                thisParentDir = obj.paramsToBackgroundDirName(thisParams);
                
            case 'Instance'
                thisParentDir = obj.paramsToInstanceDirName(thisParams);
   
            case 'threshold'
                thisParentDir = obj.paramsToThresholdDirName(thisParams);
                
            case 'psychoEllipsoid'
                thisParentDir = obj.paramsToPsychoEllipsoidDirName(thisParams);
                
            case {'Spatial'}
                thisParentDir = obj.paramsToSpatialDirName(thisParams);
                
            case {'Temporal'}
                thisParentDir = obj.paramsToTemporalDirName(thisParams);
                
            case {'Optics'}
                thisParentDir = obj.paramsToOiDirName(thisParams);
                
            case {'Mosaic'}
                thisParentDir = obj.paramsToMosaicDirName(thisParams);
            
            otherwise
                error('Unkown parent parameters type: ''%s''.', thisParams.type);
        end
        theParentDir = fullfile(theParentDir,thisParentDir);
        if (~exist(theParentDir,'dir') & p.Results.MakeDirectories)
            mkdir(theParentDir);
        else
            %fprintf(2,'Parent directory ''%s'' exists already\n', theParentDir);
        end
    end    
end

%% Filetype specific subdirectory
switch (p.Results.Type)
    case 'mat'
        theTypeDir = fullfile(theParentDir,'matfiles');
    case {'figure', 'NicePlotExport'}
        theTypeDir = fullfile(theParentDir,'figures');
    case {'movie','movieFile'}
        theTypeDir = fullfile(theParentDir,'movies');
end
if (~exist(theTypeDir,'dir') & p.Results.MakeDirectories)
    mkdir(theTypeDir);
else
    %fprintf(2,'Type directory ''%s'' exists already\n', theTypeDir);
end

%% Add in the name of reading or writing program
filedir = fullfile(theTypeDir,theProgram);
if (~exist(filedir,'dir') & p.Results.MakeDirectories)
    mkdir(filedir);
else
    %fprintf(2,'Program directory ''%s'' exists already\n', filedir);
end

%% Add in the filename
switch (p.Results.Type)
    case 'mat'
        filename = [name '.mat'];
    case {'figure', 'NicePlotExport'}
        filename = [name '.' p.Results.FigureType];
    case {'movie','movieFile'}
        filename = [name '.' p.Results.MovieType];
end
fileid = fullfile(filedir,filename);

end