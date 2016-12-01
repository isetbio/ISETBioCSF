function write(obj,name,data,paramsList,theProgram,varargin)
% write(obj,name,data,paramsList,theProgram,varargin)
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
p.addRequired('paramsList',@iscell);
p.addRequired('theProgram',@ischar);
p.addParameter('Type','mat',@ischar);
p.addParameter('ArtifactParams',[],@isstruct);
p.addParameter('FigureType','pdf',@ischar);
p.addParameter('FigureHandle', []);
p.addParameter('MovieType','m4v',@ischar);
p.parse(name,data,paramsList,theProgram,varargin{:});

% Sometimes, for compatibility, we don't actually have anything
% to write.  Handle that case.
if (isempty(data))
    return;
end

%% Get fileid
[fileid,filedir,filename] = obj.getid(p.Results.name,p.Results.paramsList,p.Results.theProgram,varargin{:},'MakeDirectories',true);

warndlg(sprintf('FileDir ''%s''.', filedir), sprintf('Filename: ''%s''', filename));

%% Write the data
switch (p.Results.Type)
    case 'mat'
        theData = p.Results.data;
        save(fileid,'theData','-v7.3');
    case 'figure'
        % The cd method seems to prevent an error when fileid gets very long.
        curdir = pwd;
        cd(filedir);
        if (exist('FigureSave','file'))
            FigureSave(filename,p.Results.data,p.Results.FigureType);
        else
            saveas(data,filename,p.Results.FigureType);
        end
        cd(curdir);
    case 'NicePlotExport'
        % The cd method seems to prevent an error when fileid gets very long.
        curdir = pwd;
        cd(filedir);
        NicePlot.exportFigToPNG(filename, p.Results.FigureHandle, 300);
        fprintf('Figure exported to %s/%s\n', pwd,filename);
        cd(curdir);
    case 'movieFile'
        % This cp followed by em seems to prevent a permissions error when the dropbox
        % location is pointed to via a softlink
        unix(['cp ' p.Results.data ' ' fileid]);
        unix(['rm ' p.Results.data]);
end

