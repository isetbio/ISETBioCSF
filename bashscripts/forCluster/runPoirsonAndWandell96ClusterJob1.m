disp('Hello !');
addpath(genpath('/home1/c/cottaris/documents/MATLAB/ToolboxToolbox'));
tbUseProject('IBIOColorDetect');
configID = 1;
c_PoirsonAndWandell96RunSession(configID, 1024);
exit;
