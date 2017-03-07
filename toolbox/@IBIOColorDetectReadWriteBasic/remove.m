function remove(obj,name,paramsList,theProgram,varargin)
% remove(obj,name,paramsList,theProgram,varargin)
%
% Delete a data object.
%
% Key/value pairs
%   'Type'
%     'mat' - Matlab .mat file
%     'figure' - A figure
%     'movie' - A movie
%     'movieFile' - String with full path to temp movie file

%% Parse input
p = inputParser;
p.addRequired('name',@ischar);
p.addRequired('paramsList',@iscell);
p.addRequired('theProgram',@ischar);
p.addParameter('Type','mat',@ischar);
p.parse(name,paramsList,theProgram,varargin{:});

%% Delete the object
fileid = obj.getid(p.Results.name,p.Results.paramsList,p.Results.theProgram,varargin{:},'MakeDirectories',false);
unix(['rm ' fileid]);

end