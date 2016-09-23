function dirname = paramsToGaborDirName(obj,spotParams)
% dirname = paramsToGaborDirName(obj,spotParams)
% 
% Generate a directory names that captures the basic mixture of monochromatic
% light spot stimulus parameters.
%
% Currently does not include background wavelength parameters

if (~strcmp(spotParams.type,'Spot'))
    error('Incorrect parameter type passed');
end

dirname = sprintf('ssz%0.0f_bsz%0.1f_swl%%d_sirrad%0.2g',...
    spotParams.spotSizeDegs,...
    spotParams.backgroundSizeDegs,...
    spotParams.spotWavelengthNm,...
    spotParams.spotCornealIrradianceUW ...
    );