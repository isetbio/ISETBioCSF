function bipolarParams = bipolarParamsGenerate()
% bipolarParams = bipolarParamsCreate(varargin)

bipolarParams.cellType = 'onDiffuse';
% sets filter as theoretical, mean physiology, or individual phys:
bipolarParams.filterType = 1; 

% sets linear, on half-wave rectification, or on and off half-wave rect
bipolarParams.rectifyType = 1;

% bp.bipolarSet('sRFcenter',1);
% bp.bipolarSet('sRFsurround',1);

% bp.bipolarSet('sRFcenter',[0 0 0; 0 1 0; 0 0 0]);
% bp.bipolarSet('sRFsurround',[0 0 0; 0 1 0; 0 0 0]);

