function [spectrum2D, sfSupport, tfSupport] = xtSpectrumFromXYTspectrum(spd3D, sfSupport, tfSupport)

    spatiotemporalSpectrumMethod = 'select sfY=0';
    spatiotemporalSpectrumMethod = 'sum over sfY';
    spatiotemporalSpectrumMethod = 'radial averaging';
    
    % organization of spd3D is [TF, SFy, SFx]
    
    switch (spatiotemporalSpectrumMethod)
        case 'sum over sfY'
            % Sum over the ySF axis
            xtSPD = squeeze(sum(spd3D,2));

            % Average along 4 quadrants
            [spectrum2D, sfSupport, tfSupport] = ...
                averageAcrossQuadrantsSpectum(sfSupport, tfSupport, xtSPD);

        case 'select sfY=0'
            % Select the ySF = 0 slice
            targetSF = 0;
            [~,sfYIndex] = min(abs(sfSupport-targetSF));
            xtSPD = squeeze(spd3D(:,sfYIndex,:));
            
            % Average along 4 quadrants
            [spectrum2D, sfSupport, tfSupport] = ...
                averageAcrossQuadrantsSpectum(sfSupport, tfSupport, xtSPD);
    
        case 'radial averaging' 
            m = numel(sfSupport)/2+1;
            sfSupport = sfSupport(m:end);
            m = numel(tfSupport)/2+1;
            tfSupport = tfSupport(m:end);
    
            indexOfZeroTF = floor(size(spd3D,1)/2)+1;
            tfIndices = indexOfZeroTF:size(spd3D,1);

            for tf = 1:numel(tfIndices)
                tfIndex = tfIndices(tf);
                tfIndexSymmetric = indexOfZeroTF - (tfIndex - indexOfZeroTF);
                if (tf == 1)
                    spectrumXY = squeeze(spd3D(tfIndex,:,:));
                else
                    spectrumXY = squeeze(spd3D(tfIndex,:,:)) + ...
                             squeeze(spd3D(tfIndexSymmetric,:,:));
                end
                meanSlice = radiallyAveragedSpectrum(spectrumXY);
                
                if (tf == 1)
                    spectrum2D = zeros(numel(tfIndices), length(meanSlice));
                end
                spectrum2D(tf,:) = meanSlice;
            end
        otherwise
            error('Unknown XT spectral method: ''%s''.', spatiotemporalSpectrumMethod)
    end
end
