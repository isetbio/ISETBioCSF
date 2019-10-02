function write(obj, name, data, paramsList, theProgram, varargin)
% Write data objects, including metadata.
%
% Syntax:
%   write(obj, name, data, paramsList, theProgram, [varargin])
%
% Description:
%    Write data objects in a good place, with metadata.
%
% Inputs:
%    obj            - Object. The IBIOColorDetectReadWriteBasic object.
%    name           - String. The object's name.
%    data           - Struct. The data structure.
%    paramsList     - Cell. A cell array of parameters and their values.
%    theProgram     - String. The function name that called write.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    Type           - String. The object type. Options are:
%           mat: Matlab .mat file. Default.
%           figure: A figure
%           movieFile: String with full path to temp movie file
%           NicePlotExportPNG: Export the plot as a PNG
%           NicePlotExportPDF: Export the plot as a PDF
%    ArtifactParams - Struct. A structure containing information specific
%                     to an artifact.
%    FigureType     - String. The type of figure to save, as supported by
%                     FigureSave/saveas. Options are:
%           pdf: PDF. Default.
%    MovieType      - String. The extension for movie filename. Options:
%           m4v: MPEG-4. Default.
%    ShowFilePath   - Boolean. Whether to display file path. Default false.
%    FigureHandle   - Handle. The figure handle. Default [].
%

%% Parse input
p = inputParser;
p.addRequired('name', @ischar);
p.addRequired('data');
p.addRequired('paramsList', @iscell);
p.addRequired('theProgram', @ischar);
p.addParameter('Type', 'mat', @ischar);
p.addParameter('ArtifactParams', [], @isstruct);
p.addParameter('FigureType', 'pdf', @ischar);
p.addParameter('FigureHandle', []);
p.addParameter('MovieType', 'm4v', @ischar);
p.addParameter('ShowFilePath', false, @islogical);
p.parse(name, data, paramsList, theProgram, varargin{:});

% Sometimes, for compatibility, we don't actually have anything
% to write.  Handle that case.
if (isempty(data)), return; end

%% Get fileid
[fileid, filedir, filename] = obj.getid(p.Results.name, ...
    p.Results.paramsList, p.Results.theProgram, varargin{:}, ...
    'MakeDirectories', true);

if (p.Results.ShowFilePath)
    warndlg(sprintf('FileDir ''%s''.', filedir), ...
        sprintf('Filename: ''%s''', filename));
end

%% Write the data
switch (p.Results.Type)
    case 'mat'
        theData = p.Results.data;
        % save(fileid, 'theData', '-v7.3');
        save(fileid, 'theData');
    case 'figure'
        % Note: The cd method seems to prevent an error when fileid gets
        % very long.
        curdir = pwd;
        cd(filedir);
        if (exist('FigureSave', 'file'))
            FigureSave(filename, p.Results.data, p.Results.FigureType);
        else
            saveas(data, filename, p.Results.FigureType);
        end
        cd(curdir);
    case 'NicePlotExportPNG'
        % Note: The cd method seems to prevent an error when fileid gets
        % very long.
        curdir = pwd;
        cd(filedir);
        NicePlot.exportFigToPNG(filename, p.Results.FigureHandle, 300);
        fprintf('Figure exported to %s/%s\n', pwd, filename);
        cd(curdir);
    case 'NicePlotExportPDF'
        % Note: The cd method seems to prevent an error when fileid gets
        % very long.
        curdir = pwd;
        cd(filedir);
        NicePlot.exportFigToPDF(filename, p.Results.FigureHandle, 300);
        fprintf('Figure exported to %s/%s\n', pwd, filename);
        cd(curdir);
    case {'movieFile', 'movie'}
        % This cp followed by em seems to prevent a permissions error when
        % the dropbox location is pointed to via a softlink
        unix(['cp ' p.Results.data ' ' fileid]);
        unix(['rm ' p.Results.data]);
end

end
