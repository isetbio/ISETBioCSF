function [theEMpaths, theEMpathsMicrons] = colorDetectMultiTrialEMPathGenerate(theConeMosaic, nTrials, eyeMovementsPerTrial, emPathType, varargin)
% theEMpaths = colorDetectEMPathGenerate(theConeMosaic, nTrials, eyeMovementsPerTrial, emPathType, varargin)
%
% Create an array of EMpaths, one for each of the nTrials, depending on the emPathType parameter
%
%  11/30/16  npc  Wrote it.
%  3/39/18   npc  Update it to use with the @fixationalEM class.
%

p = inputParser;
p.addParameter('seed',1, @isnumeric);  
p.addParameter('centeredEMPaths',false, @islogical);   
p.parse(varargin{:});

% Instantiate a fixational eye movement object for generating
% fixational eye movements that include drift and microsaccades.
fixEMobj = fixationalEM();
fixEMobj.setDefaultParams();
    
switch (emPathType)
    case 'none'
        theEMpaths = zeros(nTrials, eyeMovementsPerTrial, 2); 
        theEMpathsMicrons = theEMpaths;
    case {'frozen', 'frozen0'}
        rng(p.Results.seed);
        
        fixEMobj.computeForConeMosaic(theConeMosaic, eyeMovementsPerTrial, ...
            'nTrials', 1, ...
            'rSeed', p.Results.seed, ...
            'useParfor', false);
        theFixedEMpath = fixEMobj.emPos;
        theFixedEMpathMicrons = fixEMobj.emPosMicrons;
        
        if (strcmp(emPathType, 'frozen0'))
            k = 0;
        else
            k = 1;
        end
        theEMpaths = k*repmat(theFixedEMpath, [nTrials 1 1]);
        theEMpathsMicrons = k*repmat(theFixedEMpathMicrons, [nTrials 1 1]);
        
    case 'random'
        fixEMobj.computeForConeMosaic(theConeMosaic, eyeMovementsPerTrial, ...
            'nTrials', nTrials, ...
            'rSeed', p.Results.seed, ...
            'useParfor', false);
        theEMpaths = fixEMobj.emPos;
        theEMpathsMicrons = fixEMobj.emPosMicrons;
        
    case 'randomNoSaccades'
        fixEMobj.microSaccadeType = 'none';
        fixEMobj.computeForConeMosaic(theConeMosaic, eyeMovementsPerTrial, ...
            'nTrials', nTrials, ...
            'rSeed', p.Results.seed, ...
            'useParfor', false);
        theEMpaths = fixEMobj.emPos;
        theEMpathsMicrons = fixEMobj.emPosMicrons;
        
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
