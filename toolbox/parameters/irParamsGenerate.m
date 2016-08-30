function irParams = irParamsGenerate(varargin)
%  params.irParams = irParamsGenerate(varargin{:});
% 
% Define the parameters for the inner retina object within the IBIO
% ColorDetect tutorial framework.
% 
%   name - name of the particular IR instantiation
%   eye side - left or right
%   eye radius - radius of retinal patch in mm
%   eye ange - polar angle in degrees, 0 deg = 12 o'clock

irParams.name      = 'Macaque inner retina 1'; % This instance
irParams.eyeSide   = 'left';   % Which eye
irParams.eyeRadius = 1;        % Radius in mm
irParams.eyeAngle  = 90;       % Polar angle in degrees

irParams.runFlag = true;