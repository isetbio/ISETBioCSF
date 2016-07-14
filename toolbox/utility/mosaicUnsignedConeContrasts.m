function [LMSContrasts] = mosaicUnsignedConeContrasts(mosaicResponse,theMosaic)
% [LMSContrasts] = mosaicUnsignedConeContrasts(mosaicResponse,theMosaic)
% 
% Find the unsigned cone contrasts, given the mosaicResponses and the
% mosaic object.

% Get min max for LMS cone absorptions
% Extract the min and max absorptions in a loop. Since we are
% extracting only L, M, or S absorptions at each iteration, we get a vector
% so one call to max/min will suffice.
mosaicPattern = theMosaic.pattern;
LMSContrasts = zeros(3,1);
for ii = 2:4
    maxAbsorption = max(mosaicResponse(mosaicPattern==ii));
    minAbsorption = min(mosaicResponse(mosaicPattern==ii));
    LMSContrasts(ii-1) = (maxAbsorption-minAbsorption)/(maxAbsorption+minAbsorption);
end
