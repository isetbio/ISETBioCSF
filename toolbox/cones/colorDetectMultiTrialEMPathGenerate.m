function [theEMpaths, theEMpathsMicrons] = colorDetectMultiTrialEMPathGenerate(theConeMosaic, nTrials, eyeMovementsPerTrial, emPathType, varargin)
% theEMpaths = colorDetectEMPathGenerate(theConeMosaic, nTrials, eyeMovementsPerTrial, emPathType, varargin)
%
% Create an array of EMpaths, one for each of the nTrials, depending on the emPathType parameter
%
%  11/30/16  npc Wrote it.
%

p = inputParser;
p.addParameter('seed',1, @isnumeric);  
p.addParameter('centeredEMPaths',false, @islogical);   
p.parse(varargin{:});

switch (emPathType)
    case 'none'
        theEMpaths = zeros(nTrials, eyeMovementsPerTrial, 2); 
        theEMpathsMicrons = theEMpaths;
    case 'frozen'
        rng(p.Results.seed);
        [theFixedEMpath, theFixedEMpathMicrons] = theConeMosaic.emGenSequence(eyeMovementsPerTrial);
        theEMpaths = repmat(theFixedEMpath, [nTrials 1 1]);
        theEMpathsMicrons = repmat(theFixedEMpathMicrons, [nTrials 1 1]);
    case 'frozen0'
        rng(p.Results.seed);
        [theFixedEMpath, theFixedEMpathMicrons] = theConeMosaic.emGenSequence(eyeMovementsPerTrial);
        theEMpaths = 0*repmat(theFixedEMpath, [nTrials 1 1]);
        theEMpathsMicrons = 0*repmat(theFixedEMpathMicrons, [nTrials 1 1]);
    case 'random'
        theEMpaths = zeros(nTrials, eyeMovementsPerTrial, 2);
        for iTrial= 1:nTrials
            [theEMpaths(iTrial, :,:), theEMpathsMicrons(iTrial, :,:)] = theConeMosaic.emGenSequence(eyeMovementsPerTrial);
        end
    otherwise
        error('Unknown emPathType: ''%s''. Valid choices: ''none'', ''frozen'', ''random''.', emPathType);
end

if (p.Results.centeredEMPaths)
    % find centers
    emPathCenters = mean(theEMpaths,2);
    % subtract them from the emPaths
    theEMpaths = bsxfun(@minus, theEMpaths, emPathCenters);
    % emPaths must be integered-valued, so round
    theEMpaths = round(theEMpaths);
    
    % find centers
    emPathCenters = mean(theEMpathsMicrons,2);
    % subtract them from the emPaths
    theEMpathsMicrons = bsxfun(@minus, theEMpathsMicrons, emPathCenters);    
end
end
