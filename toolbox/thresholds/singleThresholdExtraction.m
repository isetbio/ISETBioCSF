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
options.MaxFunEvals = 1000;
options.MaxIter     = 1000;
options.Display     = 'off';

%% Search grid specification for Palemedes
gridLevels = 100;
searchGrid.alpha = logspace(log10(stimLevels(1)),log10(stimLevels(end)),gridLevels);
searchGrid.beta = 10.^linspace(-4,4,gridLevels);
searchGrid.gamma = 0.5;
searchGrid.lambda = 0.0;

%% Use Palamedes grid search method
[paramsValues,LL,flag] = PAL_PFML_Fit(stimLevels(:), numPos(:), outOfNum(:), ...
            searchGrid, paramsFree, PF, 'SearchOptions', options);

%% Get threshold and deal with catastrophic cases
threshold = PF(paramsValues, criterionFraction, 'inverse');
if (threshold < 0 | threshold > 1 | ~isreal(threshold) | isinf(threshold))
    threshold = NaN;
end

%% Provide fit psychometric function on passed stimulus levels
if (~isnan(threshold))
    fitFractionCorrect = PF(paramsValues,fitStimLevels);
else
    fitFractionCorrect = NaN*ones(size(fitStimLevels));
end

end

