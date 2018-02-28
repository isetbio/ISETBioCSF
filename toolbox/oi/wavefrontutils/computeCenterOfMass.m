function centerOfMass = computeCenterOfMass(PSF, varargin)
    p = inputParser;
    p.addParameter('useMaxAsCenterOfMass', true, @islogical);
    p.parse(varargin{:});
    
    if (p.Results.useMaxAsCenterOfMass)
        [~,idx] = max(PSF(:));
        [i,j] = ind2sub(size(PSF), idx);
        centerOfMass = [j i];
    else
        [rc,cc] = ndgrid(1:size(PSF,1),1:size(PSF,2));
        Mt = sum(PSF(:));
        centerOfMassY = sum(PSF(:) .* rc(:)) / Mt;
        centerOfMassX = sum(PSF(:) .* cc(:)) / Mt;
        centerOfMass = [centerOfMassX centerOfMassY];
    end
end
