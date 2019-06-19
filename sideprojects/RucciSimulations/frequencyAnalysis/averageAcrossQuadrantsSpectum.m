function [averageSpectra, sfAxis, tfAxis] = averageAcrossQuadrantsSpectum(sfSupport, tfSupport, xtSpectra)

    indexOfZeroSF = floor(size(xtSpectra,2)/2)+1;
    indexOfZeroTF = floor(size(xtSpectra,1)/2)+1;
    
    sfIndices = indexOfZeroSF:size(xtSpectra,2);
    tfIndices = indexOfZeroTF:size(xtSpectra,1);
    averageSpectra = zeros(numel(tfIndices), numel(sfIndices));
    
    sfAxis = sfSupport(sfIndices);
    tfAxis = tfSupport(tfIndices);

    for tf = 1:numel(tfIndices)
        tfIndex = tfIndices(tf);
        tfIndexSymmetric = indexOfZeroTF - (tfIndex - indexOfZeroTF);
        for sf = 1:numel(sfIndices)
            sfIndex = sfIndices(sf);
            sfIndexSymmetric = indexOfZeroSF - (sfIndex - indexOfZeroSF);
            if ((tf == 1) && (sf == 1))
                % (0,0) frequency
                averageSpectra(tf,sf) = xtSpectra(tfIndex ,sfIndex);
            else
                if (tf == 1) && (sf > 1)
                    % at zero TF axis, sum corresponding points from 2
                    % spatial frequencies
                    averageSpectra(tf,sf) = ...
                        xtSpectra(tfIndex,sfIndex) + ...
                        xtSpectra(tfIndex,sfIndexSymmetric);
                elseif (sf == 1) && (tf > 1)
                    % at zero SF axis, sum corresponding points from 2
                    % temporal frequencies
                    averageSpectra(tf,sf) = ...
                        xtSpectra(tfIndex,sfIndex) + ...
                        xtSpectra(tfIndexSymmetric,sfIndex);
                else
                    % sum corresponding points from 4 quadrants
                    averageSpectra(tf,sf) = ...
                        xtSpectra(tfIndex,sfIndex) + ...
                        xtSpectra(tfIndex,sfIndexSymmetric) + ...
                        xtSpectra(tfIndexSymmetric,sfIndex) + ...
                        xtSpectra(tfIndexSymmetric,sfIndexSymmetric);
                end
            end
        end
    end 
end
