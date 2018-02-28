%availableCustomWvfOpticsModels   Return cell array with names of available custom wvf optic
%   availableCustomWvfOpticsModels()
%
% Also see: oiWithCustomOptics
%
% 6/20/17  npc  Wrote it.

function opticsModels = availableCustomWvfOpticsModels

    % Available custom wvf optics models
    opticsModels = {...
        'AOoptics75mmPupil' ...
        'AOoptics80mmPupil' ...
        'WvfHuman' ...
        'WvfHumanMeanOTFmagMeanOTFphase' ...        % Mean OTF mag and phase (over 100 subjects)
        'WvfHumanMeanOTFmagZeroOTFphase' ...        % Mean OTF mag, zero OTF phase(over 100 subjects)
        'ThibosBestPSFSubject3MMPupil' ...          %  
        'ThibosDefaultSubject3MMPupil' ...          %  
        'ThibosAverageSubject3MMPupil' ...          %  
        'ThibosDefocusedSubject3MMPupil', ...       %  
        'ThibosVeryDefocusedSubject3MMPupil' ...    %  
    };

end

