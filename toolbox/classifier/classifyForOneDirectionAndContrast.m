function [usePercentCorrect,useStdErr,h] = classifyForOneDirectionAndContrast(noStimData,stimData,thresholdParams,varargin)
% [usePercentCorrect,useStdErr,h] = classifyForOneDirectionAndContrast(noStimData,stimData,thresholdParams,varargin)
%
% Do classification for one stimulus color direction and contrast.
%
% Optional key/value pairs
%   'plotSvmBoundary' - true/false (default false).  Plot classification boundary
%   'plotPCAAxis1' - First PCA component to plot (default 1)
%   'plotPCAAxis2' - Second PCA component to plot (default 2)
%
% 7/14/16  dhb, xd  Pulled this out.  Hoping it will work.

% Parse args
p = inputParser;
p.addRequired('noStimData',@isstruct);
p.addRequired('stimData',@isstruct);
p.addRequired('thresholdParams',@isstruct);
p.addParameter('plotSvmBoundary',false,@islogical);
p.addParameter('plotPCAAxis1',1,@isnumeric)
p.addParameter('plotPCAAxis2',2,@isnumeric)
p.parse(noStimData,stimData,thresholdParams,varargin{:});

%% Put zero contrast response instances into data that we will pass to the SVM
[classificationData,classes] = classificationDataNoStimDataInitialize(noStimData,thresholdParams);

%% Insert stimulus responses into data for classifier
[classificationData,classes] = classificationDataStimDataInsert(classificationData,classes,stimData,thresholdParams);

%% Decide what type of classifier we are running
switch (thresholdParams.method)
    case 'svm'
        % Friendly neighborhood SVM, with optional standardization and PCA
        % first
        
        % Do PCA?  This also standardizes the data if desired, and gets rid of any
        % features with no variance.
        fprintf('\tDoing PCA etc ... ');
        theData = transformDataWithPCA(classificationData,thresholdParams.PCAComponents,thresholdParams.STANDARDIZE);
        fprintf('done\n');
        
        % Perform SVM classification for this stimulus vs the zero contrast stimulus
        fprintf('\tRunning SVM ...');
        [usePercentCorrect, useStdErr, svm] = classifyWithSVM(theData,classes,thresholdParams.kFold);
        fprintf(' correct: %2.2f%%\n', usePercentCorrect*100);
        
        % Optional plot of classification boundary
        if (p.Results.plotSvmBoundary)
            % Plot the data points along these two components
            PCAAxis1 = p.Results.plotPCAAxis1;
            PCAAxis2 = p.Results.plotPCAAxis2;
            
            h = figure('Position',[0 0 750 750]);
            hold on;
            indexA = find(classes == 0);
            indexB = find(classes == 1);
            plot(theData(indexA,PCAAxis1),theData(indexA,PCAAxis2),'*','MarkerSize',7);
            plot(theData(indexB,PCAAxis1),theData(indexB,PCAAxis2),'o','MarkerSize',7);
            
            % Extract the boundary line. svm.Beta is the line orthogonal to the
            % boundary. We take it's null space (and since we only have 2 dimensions,
            % it is a single vector) which is the boundary. We scale it by a
            % large number (100) so the line comes out nicely.
            orthLineToBoundary = svm.Beta;
            orthLineToBoundary = orthLineToBoundary([PCAAxis1 PCAAxis2]);
            boundary = null(orthLineToBoundary')*100;
            
            % Plot the decision boundary
            plot([-boundary(1) boundary(1)],[-boundary(2) boundary(2)],'--k','LineWidth',2);
            
            % Reset the axis limits
            ylim([min(theData(:,PCAAxis2)) max(theData(:,PCAAxis2))]);
            xlim([min(theData(:,PCAAxis1)) max(theData(:,PCAAxis1))]);
            
            % Make plot look nice + various labels/titles/legend
            set(gca,'FontName','Helvetica','FontSize',18,'LineWidth',2);
            legend({'Null Stimulus','Stimulus'},'Location','Northeast','FontSize',18);
            
            axis square;
            xlabel(sprintf('Principal Component %d',PCAAxis1),'FontSize',22);
            ylabel(sprintf('Principal Component %d',PCAAxis2),'FontSize',22);
            title('SVM Decision Boundary on Two Principal Components','FontSize',24);
        else
            h = [];
        end
           
    case 'mlpt'
        % Template maximum likelihood classifier, based on analytic mean responses o
        % each class and assuming that the noise is Poisson.  
        fprintf('\tRunning template LogLikelihood classifier ...\n');
        
        % Check for sanity
        if (~strcmp(thresholdParams.signalSource,'isomerizations'))
            error('Template only available at isomerizations');
        end
     
        % Set up template
        if (thresholdParams.nIntervals == 1)
            templateClass0 = noStimData.noiseFreeIsomerizations(:)';
            templateClass1 = stimData.noiseFreeIsomerizations(:)';
        else
            templateClass0 = [noStimData.noiseFreeIsomerizations(:)' stimData.noiseFreeIsomerizations(:)'];
            templateClass1 = [stimData.noiseFreeIsomerizations(:)' noStimData.noiseFreeIsomerizations(:)'];
        end
        
        % Compute likelihood of each class, given empirical means and use
        % this to predict class
        teClasses = classes;
        teData = classificationData;
        tePredict = -1*ones(size(teClasses));
        nTeObservations = size(teData,1);
        for ii = 1:nTeObservations
            loglGiven0(ii) = sum(log10(poisspdf(teData(ii,:),templateClass0)));
            loglGiven1(ii) = sum(log10(poisspdf(teData(ii,:),templateClass1)));
            if (loglGiven1(ii) > loglGiven0(ii))
                tePredict(ii) = 1;
            else
                tePredict(ii) = 0;
            end
        end
        
        nCorrect = length(find(tePredict == teClasses));
        usePercentCorrect = nCorrect/length(teClasses);
        useStdErr = 0;
        
        % Report percent correct
        fprintf('\tPercent correct: %2.2f%%\n\n', usePercentCorrect*100);
        
        % No figure with this method
        h = [];  
        
    case 'mlpe'
        % Empirical Poisson maximum likelihood classifier, based on mean responses to
        % each class and assuming that the noise is Poisson.
        fprintf('\tRunning empirical LogLikelihood classifier ...\n');
        
        % Create cross-validation partition
        nObservations = size(classificationData,1);
        crossVal = cvpartition(nObservations,'KFold',thresholdParams.kFold);
        
        % Loop through cross validations
        nCrossVals = crossVal.NumTestSets;
        for cc = 1:nCrossVals
            fprintf('\t\tCross validation %d of %d\n',cc,nCrossVals);
            trIdx = crossVal.training(cc);
            teIdx = crossVal.test(cc);
            
            % Get mean response of each class from training set
            trClasses = classes(trIdx);
            trData = classificationData(trIdx,:);
            
            class0Index = trClasses == 0;
            class1Index = trClasses == 1;
            meanTrClass0 = mean(trData(class0Index,:),1);
            meanTrClass1 = mean(trData(class1Index,:),1);
            
            % Classify test set
            %
            % % Compute likelihood of each class, given empirical means and use
            % this to predict class
            teClasses = classes(teIdx);
            teData = classificationData(teIdx,:);
            tePredict = -1*ones(size(teClasses));
            nTeObservations = size(teData,1);
            for ii = 1:nTeObservations
                loglGiven0(ii) = sum(log10(poisspdf(teData(ii,:),meanTrClass0)));
                loglGiven1(ii) = sum(log10(poisspdf(teData(ii,:),meanTrClass1)));
                if (loglGiven1(ii) > loglGiven0(ii))
                    tePredict(ii) = 1;
                else
                    tePredict(ii) = 0;
                end
            end
        
            nCorrect = length(find(tePredict == teClasses));
            percentCorrect(cc) = nCorrect/length(teClasses);
            fprintf('\t\tPercent correct = %2.2f%%\n',100*percentCorrect(cc));
        end
        
        % Aggregate percent correct
        usePercentCorrect = mean(percentCorrect);
        useStdErr = std(percentCorrect)/sqrt(length(percentCorrect));;

        % Report percent correct
        fprintf('\tPercent correct: %2.2f%%\n\n', usePercentCorrect*100);
        
        % No figure with this method
        h = [];  
        
    otherwise
        error('Unknown classification method');
end

end