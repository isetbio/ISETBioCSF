function [percentCorrect,stdErr,svm] = classifyWithSVM(data,classes,kFold)
% [percentCorrect,stdErr,svm] = classifyWithSVM(data,classes,kFold)
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
% p = inputParser;
% p.addRequired('data',@isnumeric);
% validateClasses = @(X) numel(unique(X))==2;
% p.addRequired('classes',validateClasses);
% p.parse(data,classes);

%% Train cross validated SVM
svm = fitcsvm(data,classes,'KernelFunction','linear','KernelScale','auto');
CVSVM = crossval(svm,'KFold',kFold);
percentCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
stdErr = std(percentCorrect)/sqrt(kFold);
percentCorrect = mean(percentCorrect);

% We'll discard the data that's normally stored in the SVM object to save space.
svm = compact(svm);
svm = discardSupportVectors(svm);

end

