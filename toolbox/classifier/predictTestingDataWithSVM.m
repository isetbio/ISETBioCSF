function percentCorrect = predictTestingDataWithSVM(svm,testingData,classes)
% percentCorrect = predictTestingDataWithSVM(svm,testingData,classes)
% 
% Returns percent correct of classification of testingData by svm.
%
% 7/12/16  xd  wrote it

percentCorrect = sum(predict(svm,testingData)==classes)/length(classes);

end

