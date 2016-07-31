function [data,artifactParams] = read(obj,name,parentParamsList,currentParamsList,theProgram,varargin)
% [data,artifactParams] = read(obj,name,parentParamsList,currentParamsList,theProgram,varargin)
%
% Read data objects from the right place.
%
% Key/value pairs
%   'Type'
%     'mat' - Matlab .mat file

%% Parse input
p = inputParser;
p.addRequired('name',@ischar);
p.addRequired('parentParamsList',@iscell);
p.addRequired('currentParamsList',@iscell);
p.addRequired('theProgram',@ischar);
p.addParameter('Type','mat',@ischar);
p.parse(name,parentParamsList,currentParamsList,theProgram,varargin{:});

%% Get fileid
fileid = obj.getid(p.Results.name,p.Results.parentParamsList,p.Results.currentParamsList,p.Results.theProgram,varargin{:});

%% Read the data
switch (p.Results.Type)
    case 'mat'
        dataStruct = load(fileid,'theData');
        data = dataStruct.theData;
        artifactParams = [];
end

end