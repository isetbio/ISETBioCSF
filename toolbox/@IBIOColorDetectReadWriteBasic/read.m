function data = read(obj,fileid,varagin)
% data = read(obj,fileid,varagin)
%
% Read data objects, given unique identifier
%
% Key/value pairs
%   'Type'
%     'mat' - Matlab .mat file (default)
%     'figure' - A figure
%     'movie' - A movie
%     'movieFile' - String with full path to temp movie file

% Parse input.
p = inputParser;
p.addRequired('name',@ischar);
p.addParameter('Type','mat',@ischar);

p.parse(params,fileid,varargin{:});
fileid = p.Results.name;

end