function dirname = paramsToSpotSpatialDirName(obj,spotParams)
% dirname = paramsToSpotSpatialDirName(obj,spotParams)
% 
% Generate a directory names that captures the basic mixture of monochromatic
% light spot stimulus parameters.
%
% Currently does not include background wavelength parameters

if (~strcmp(spotParams.type,'SpotSpatial'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('ssz%0.0f_bsz%0.1f',...
    spotParams.spotSizeDegs,...
    spotParams.backgroundSizeDegs,...
    );