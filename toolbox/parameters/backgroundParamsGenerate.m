function backgroundParams = backgroundParamsGenerate(varargin)
% backgroundParams = backgroundParamsGenerate(varargin)
%
% Generate parameters for a color modulation.
%
%   backgroundxYY - Colorimetric specification of background, in CIE xyY (cd/m2)
%   monitorFile - Isetbio display description of monitor on which grating is shown.
%   leakageLum - Luminance when monitor input is zero.
%   lumFactor - Multiply background luminance by this factor.
%               Adjust display to make it work.

backgroundParams.type = 'Background';

backgroundParams.backgroundxyY = [0.27 0.30 49.8]';
backgroundParams.monitorFile = 'CRT-MODEL';
backgroundParams.leakageLum = 1.0;
backgroundParams.lumFactor = 1;