function fileid = read(obj,name,varagin)
% fileid = read(obj,name,varagin)
%
% Read data objects using metadata.
%
% Key/value pairs
%   'Type'
%     'mat' - Matlab .mat file (default)
%     'figure' - A figure
%     'movie' - A movie
%     'movieFile' - String with full path to temp movie file
%   'ArtifactParams'
%      Structure with artifact specific information.

% Parse input.
p = inputParser;
p.addRequired('name',@ischar);
p.addParameter('Type','mat',@ischar);
p.addParameter('ArtifactParams','mat',@isstruct);

p.parse(params,name,varargin{:});
name = p.Results.name;

end