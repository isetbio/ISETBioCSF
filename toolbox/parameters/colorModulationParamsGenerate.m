function colorModulationParams = colorModulationParamsGenerate(varargin)
% colorModulationParams = colorModulationParamsGenerate(varargin)
%
% Generate parameters for a color modulation.
%
%   contrast - Contrast specfied relative to coneContrasts.
%              Can be a vector of contrasts.
%   coneContrasts - Color direction of grating in cone contrast space
%                   Can be a 3 by N matrix of contrast directions.

colorModulationParams.type = 'ColorModulation';

colorModulationParams.contrast = 1;
colorModulationParams.coneContrasts = [0.05 -0.05 0]';
