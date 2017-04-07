function [percentCorrect,stdErr,svm] = classifyWithSVM(data,classes,kFold, varargin)
% [percentCorrect,stdErr,svm] = classifyWithSVM(data,classes,kFold, varargin)
% 
% Trains a SVM using kFold Cross Validation and returns the average
% percent correct. The data will be divided into ten roughly evenly sized
% sets. For each set, the SVM will be trained on the other nine and tested
% on the set held out. Then, the average percent correct is returned.
%
% Inputs:
%   data     -   A matrix containing entries of data along the rows and
%                features along the columns.
%   
%   classes  -   A vector containing class assignments. Each entry
%                corresponds to a row in the data matrix. 
%
% 7/7/16  xd  wrote it

%% Parse inputs
p = inputParser;
p.addParameter('standardizeSVMpredictors', false, @islogical);
p.addParameter('useRBFKernel', false, @islogical);
p.parse(varargin{:});
standardizeSVMpredictors = p.Results.standardizeSVMpredictors;

%% Train cross validated SVM
if (p.Results.useRBFKernel)
    kernelFunction = 'rbf';
else
    kernelFunction = 'linear';
end

svm = fitcsvm(data,classes,'KernelFunction',kernelFunction,'KernelScale','auto', 'standardize', standardizeSVMpredictors);
CVSVM = crossval(svm,'KFold',kFold);
percentCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
stdErr = std(percentCorrect)/sqrt(kFold);
percentCorrect = mean(percentCorrect);

% We'll discard the data that's normally stored in the SVM object to save space.
% This can only be done in the LINEAR SVM
svm = compact(svm);
if (~p.Results.useRBFKernel)
    svm = discardSupportVectors(svm);
end
end

