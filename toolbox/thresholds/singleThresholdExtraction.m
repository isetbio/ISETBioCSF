function [threshold,fitFractionCorrect,paramsValues] = singleThresholdExtraction(stimLevels,fractionCorrect,criterionFraction,numTrials,fitStimLevels)
% [threshold,fitFractionCorrect,paramsValues] = singleThresholdExtraction(stimLevels,fractionCorrect,criterionFraction,numTrials,fitStimLevels)
% 
% Fits a cumulative Weibull to the data variable and returns
% the threshold at the criterion as well as the parameters needed to plot the
% fitted curve. 
%
% Note that this function requires and makes use of Palamedes Toolbox 1.8
%
% Inputs:
%   stimLevels   -  A N-vector representing stimuli levels corresponding to the data.
%
%   fractionCorrect -  A N-vector of [0-1] fractions that represent performance on a task.
%
%   criterionFraction -  Fraction correct at which to find a threshold.
%
%   numTrials - Scalar, number of trials run for each stimulus level (assumed same for all levels).
%
% xd  6/21/16 wrote it

%% Set some parameters for the curve fitting
paramsFree     = [1 1 0 0];
numPos = round(numTrials*fractionCorrect);
outOfNum       = repmat(numTrials,1,length(fractionCorrect));
PF             = @PAL_Weibull;

%% Some optimization settings for the fit
options             = optimset('fminsearch');
options.TolFun      = 1e-09;
options.MaxFunEvals = 10000*100;
options.MaxIter     = 500*100;
options.Display     = 'off';

%% Fit psychometric function to the data
%
% Try starting the search at various initial values and keep the best.
tryIndex = 1;
startSlopes = [0.001 0.01 0.1 1];
startLocations = stimLevels;
for jj = 1:length(startSlopes)
    for ii = 1:length(startLocations)
        paramsEstimate = [startLocations(ii) startSlopes(jj) 0.5 0];
        [paramsValuesAll{tryIndex},LLAll(tryIndex)] = PAL_PFML_Fit(stimLevels(:), numPos(:), outOfNum(:), ...
            paramsEstimate, paramsFree, PF, 'SearchOptions', options);
        thresholdAll(tryIndex) = PF(paramsValuesAll{tryIndex}, criterionFraction, 'inverse');
        if (paramsValuesAll{tryIndex}(1) < 0 | paramsValuesAll{tryIndex}(2) < 0)
            LLAll(tryIndex) = -Inf;
        end
        tryIndex = tryIndex + 1;       
    end
end
[~,bestIndex] = max(real(LLAll));
bestII = bestIndex(1);
paramsValues = paramsValuesAll{bestII};
threshold = thresholdAll(bestII);

% Deal with catostrophic cases
if (threshold > 1 | ~isreal(threshold) | isinf(threshold))
    threshold = NaN;
end

% Provide fit psychometric function on passed stimulus levels
fitFractionCorrect = PF(paramsValues,fitStimLevels);

end

