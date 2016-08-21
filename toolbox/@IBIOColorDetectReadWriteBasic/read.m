function [data,artifactParams] = read(obj,name,paramsList,theProgram,varargin)
% [data,artifactParams] = read(obj,name,paramsList,theProgram,varargin)
%
% Read data objects from the right place.
%
% Key/value pairs
%   'Type'
%     'mat' - Matlab .mat file

%% Parse input
p = inputParser;
p.addRequired('name',@ischar);
p.addRequired('paramsList',@iscell);
p.addRequired('theProgram',@ischar);
p.addParameter('Type','mat',@ischar);
p.parse(name,paramsList,theProgram,varargin{:});

%% Get fileid
fileid = obj.getid(p.Results.name,p.Results.paramsList,p.Results.theProgram,varargin{:},'MakeDirectories',false);

%% Read the data
switch (p.Results.Type)
    case 'mat'
        dataStruct = load(fileid,'theData');
        data = dataStruct.theData;
        artifactParams = [];
end

end