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
%   'MovieType'
%      Extension for movie filename
%      'm4v' - MPEG-4 (default)

%% Get fileid
fileid = obj.getid(name,varargin{:});

%% Write the data
switch (p.Results.Type)
    case 'mat'
        save(fileid),data,'-v7.3');
    case 'figure'
        if (exist('FigureSave','file'))
            FigureSave(fullfile(fileid,data,p.Results.FigureType);
        else
            saveas(h,fileid,p.Results.FigureType);
        end
    case 'movieFile'
        unix(['mv ' data ' ' fileid]);
end

