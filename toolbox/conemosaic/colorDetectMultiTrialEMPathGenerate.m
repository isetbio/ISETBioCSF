function theEMpaths = colorDetectMultiTrialEMPathGenerate(theConeMosaic, nTrials, eyeMovementsPerTrial, emPathType, varargin)
% theEMpaths = colorDetectEMPathGenerate(theConeMosaic, nTrials, eyeMovementsPerTrial, emPathType, varargin)
%
% Create an array of EMpaths, one for each of the nTrials, depending on the emPathType parameter
%
%  11/30/16  npc Wrote it.
%

p = inputParser;
p.addParameter('seed',1, @isnumeric);   
p.parse(varargin{:});

switch (emPathType)
    case 'none'
        theEMpaths = zeros(nTrials, eyeMovementsPerTrial, 2);     
    case 'frozen'
        rng(p.Results.seed);
        theFixedEMpath = theConeMosaic.emGenSequence(eyeMovementsPerTrial);
        theEMpaths = permute(repmat(theFixedEMpath, [1 1 nTrials]), [3 1 2]);
    case 'random'
        theEMpaths = zeros(nTrials, eyeMovementsPerTrial, 2);
        for iTrial= 1:nTrials
            theEMpaths(iTrial, :,:) = theConeMosaic.emGenSequence(eyeMovementsPerTrial);
        end
    otherwise
        error('Unknown emPathType: ''%s''. Valid choices: ''none'', ''frozen'', ''random''.', emPathType);
end
end
