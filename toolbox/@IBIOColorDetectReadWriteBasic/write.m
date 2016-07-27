function write(obj,name,data,parentParamsList,currentParamsList,theProgram,varargin)
% write(obj,name,data,parentParamsList,currentParamsList,theProgram,varargin)
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
%   'MovieType'
%      Extension for movie filename
%      'm4v' - MPEG-4 (default)

%% Parse input
p = inputParser;
p.addRequired('name',@ischar);
p.addRequired('data');
p.addRequired('parentParamsList',@iscell);
p.addRequired('currentParamsList',@iscell);
p.addRequired('theProgram',@ischar);
p.addParameter('Type','mat',@ischar);
p.addParameter('ArtifactParams',[],@isstruct);
p.addParameter('FigureType','pdf',@ischar);
p.addParameter('MovieType','m4v',@ischar);
p.parse(name,data,parentParamsList,currentParamsList,theProgram,varargin{:});

%% Get fileid
fileid = obj.getid(p.Results.name,p.Results.parentParamsList,p.Results.currentParamsList,p.Results.theProgram,varargin{:});

%% Write the data
switch (p.Results.Type)
    case 'mat'
        theData = p.Results.data;
        save(fileid,'theData','-v7.3');
    case 'figure'
        if (exist('FigureSave','file'))
            FigureSave(fileid,p.Results.data,p.Results.FigureType);
        else
            saveas(h,fileid,p.Results.FigureType);
        end
    case 'movieFile'
        unix(['mv ' p.Results.data ' ' fileid]);
end

