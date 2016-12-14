function modStruct = modifyStructParams(oldStruct, varargin)
    modStruct = oldStruct;
    for k = 1:2:numel(varargin)
        if (isfield(modStruct, varargin{k}))
            modStruct.(varargin{k}) = varargin{k+1};
        else
        	error('field: ''%s'', does not exist in input struct.', varargin{k});
        end
    end

end

