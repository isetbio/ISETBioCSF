function colorModulationParams = colorModulationParamsGenerate(varargin)
% colorModulationParams = colorModulationParamsGenerate(varargin)
%
% Generate parameters for a color modulation.
%
%   contrast - Contrast specfied relative to coneContrasts.
%              Can be a vector of contrasts.
%   coneContrasts - Color direction of grating in cone contrast space
%                   Can be a 3 by N matrix of contrast directions.
%   backgroundxYY - Colorimetric specification of background, in CIE xyY (cd/m2)
%   monitorFile - Isetbio display description of monitor on which grating is shown.
%   leakageLum - Luminance when monitor input is zero.

colorModulationParams.type = 'ColorModulation';

colorModulationParams.contrast = 1;
colorModulationParams.coneContrasts = [0.05 -0.05 0]';
colorModulationParams.backgroundxyY = [0.27 0.30 49.8]';
colorModulationParams.monitorFile = 'CRT-MODEL';
colorModulationParams.leakageLum = 1.0;