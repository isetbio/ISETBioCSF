function bipolarParams = bipolarParamsGenerate()
% bipolarParams = bipolarParamsCreate(varargin)
% 
% cell type - type of bipolar, determines spatial pooling and polarity
% filter type - determines how temporal filter is built, should always be
%                   type 1.
% rectify type - 1 is linear center, surround set to zero
%                2 is rectified center, surround set to zero
%                3 is linear center, linear surround
%                4 is rectified center, rectified surround

bipolarParams.cellType = 'onDiffuse';
% sets filter as theoretical, mean physiology, or individual phys:
bipolarParams.filterType = 1; 

% sets linear, on half-wave rectification, or on and off half-wave rect
bipolarParams.rectifyType = 1;

% bp.bipolarSet('sRFcenter',1);
% bp.bipolarSet('sRFsurround',1);

% bp.bipolarSet('sRFcenter',[0 0 0; 0 1 0; 0 0 0]);
% bp.bipolarSet('sRFsurround',[0 0 0; 0 1 0; 0 0 0]);

