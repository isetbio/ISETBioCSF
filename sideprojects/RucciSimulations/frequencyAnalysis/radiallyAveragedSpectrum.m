function meanSlice = radiallyAveragedSpectrum(spectrum2D)

    if (ndims(spectrum2D) ~= 2)
        error('2D spatial spectrum expected')
    end
    midRow = floor(size(spectrum2D,1)/2);
    
    % make it odd, so we have a clear center for rotation
    if (mod(size(spectrum2D,1),2) == 0) && (mod(size(spectrum2D,2),2) == 0)
        % remove first row and col
        spectrum2D = spectrum2D(2:end, 2:end);
    else
        error('problem with size');
    end
    
    rotationAngles = 0:2:359;
    rotatedSlices = zeros(numel(rotationAngles), size(spectrum2D,2));
    parfor k = 1:numel(rotationAngles)
        s = imrotate(spectrum2D,rotationAngles(k),'bilinear','crop');
        rotatedSlices(k,:) = s(midRow,:);
    end
    meanSlice = sum(rotatedSlices,1);
    meanSlice = meanSlice(midRow:end);
    
%     figure(1)
%     subplot(1,2,1)
%     imagesc(spectrum2D(midRow:end, midRow:end));
%     
%     subplot(1,2,2);
%     plot(meanSlice, 'k-');
%     pause
end

