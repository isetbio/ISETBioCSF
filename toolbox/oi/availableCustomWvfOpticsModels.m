%availableCustomWvfOpticsModels   Return cell array with names of available custom wvf optic
%   availableCustomWvfOpticsModels()
%
% Also see: oiWithCustomOptics
%
% 6/20/17  npc  Wrote it.

function opticsModels = availableCustomWvfOpticsModels

    % Available custom wvf optics models
    opticsModels = {...
        'AOoptics75mmPupil', ...
        'AOoptics80mmPupil', ...
        'WvfHumanMeanOTFmagMeanOTFphase', ...       % Mean OTF mag and phase (over 100 subjects)
        'WvfHumanMeanOTFmagZeroOTFphase', ...       % Mean OTF mag, zero OTF phase(over 100 subjects)
        'WvfHumanSubject1', ...                     % OTF from selected subject #1
        'WvfHumanSubject2', ...                     % OTF from selected subject #2
        'WvfHumanSubject3', ...                     % OTF from selected subject #3
        'WvfHumanSubject4', ...                     % OTF from selected subject #4
        'WvfHumanSubject5' ...                      % OTF from selected subject #5
    };

end

