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
% 1/17/17  npc      Added svmV1FilterBank classifier

% Parse args        
p = inputParser;
p.addRequired('noStimData',@isstruct);
p.addRequired('stimData',@isstruct);
p.addRequired('thresholdParams',@isstruct);
p.addParameter('plotSvmBoundary',false,@islogical);
p.addParameter('plotPCAAxis1',1,@isnumeric);
p.addParameter('plotPCAAxis2',2,@isnumeric);
p.parse(noStimData,stimData,thresholdParams,varargin{:});

%% Transform the raw cone responses into V1 filter bank responses
if (strcmp(thresholdParams.method, 'svmV1FilterBank'))
    if ((~isfield(thresholdParams, 'V1filterBank')) || (isfield(thresholdParams, 'V1filterBank')) && (isempty(thresholdParams.V1filterBank)))
        error('thresholdParams must have a V1filterBank field when using the svmV1FilterBank classifier\n');
    end
    [noStimData, stimData] = transformDataWithV1FilterBank(noStimData, stimData, thresholdParams.V1filterBank,thresholdParams.STANDARDIZE);
end

%% Put zero contrast response instances into data that we will pass to the SVM
[classificationData,classes] = classificationDataNoStimDataInitialize(noStimData,thresholdParams);

%% Insert stimulus responses into data for classifier
[classificationData,classes] = classificationDataStimDataInsert(classificationData,classes,stimData,thresholdParams);

%% Decide what type of classifier we are running
switch (thresholdParams.method)
    
    case 'svmV1FilterBank'
        % Perform SVM classification for this stimulus vs the zero contrast stimulus
        fprintf('\tRunning SVM ...');
        [usePercentCorrect, useStdErr, svm] = classifyWithSVM(classificationData,classes,thresholdParams.kFold,thresholdParams.standardizeSVMpredictors);
        fprintf(' correct: %2.2f%%\n', usePercentCorrect*100);
        h = [];
        
    case 'svm'
        % Friendly neighborhood SVM, with optional standardization and PCA
        % first
        
        % Do PCA?  This also standardizes the data if desired, and gets rid of any
        % features with no variance.
        fprintf('\tExtracting the first %d principal components of the data ... ', thresholdParams.PCAComponents);
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
        
        % Originally in this code, the likelihood ratio was computed via
        % the poisspdf function.  Why not, you'd think?  But, this turns
        % out to have numerical problems, or something.  Things work
        % considerably better if you take advantage of knowing the analytic
        % form of the Poisson PDF and working out the form of the log
        % likelihood ratio.
        %
        % The problem with the poisspdf method manifested itself in lower
        % percent than predicted analytically, and this in turn was related
        % to 0 not being the best criterion.
        %
        % When you look at the form of the Poisson PDF, you realize that
        % there may well be some approximations made for large N, because
        % computing the various factorials and exponentials is probably
        % problematic.  Going to the log ratio directly eliminates the
        % middleman, and in particular the factorial cancels and goes away.
        %
        % A nice side effect is that the analytic method is considerably
        % faster.
        
        % We only need to consider data where the two classes differ, so
        % find that subset and use below.
        index = find(templateClass0 ~= templateClass1);
        ANALYTIC_LIKELY = true;
        if (ANALYTIC_LIKELY)
            logTemplateDiff = log(templateClass1(index))-log(templateClass0(index));
            C = sum(templateClass0(index)-templateClass1(index));
            for ii = 1:nTeObservations
                logLikelyRatio(ii) = sum( teData(ii,index) .* logTemplateDiff ) + C;
            end
        else
            for ii = 1:nTeObservations       
                loglGiven0(ii) = sum(log10(poisspdf(teData(ii,index),templateClass0(index))));
                loglGiven1(ii) = sum(log10(poisspdf(teData(ii,index),templateClass1(index))));
                logLikelyRatio(ii) = loglGiven1(ii)-loglGiven0(ii);
            end
        end
        
        % Compute percent correct
        for ii = 1:nTeObservations
            if (logLikelyRatio(ii) > 0)
                tePredict(ii) = 1;
            else
                tePredict(ii) = 0;
            end
        end 
        nCorrect = length(find(tePredict == teClasses));
        usePercentCorrect = nCorrect/length(teClasses);
        useStdErr = 0;
        fprintf('\tPercent correct: %2.2f%%\n',usePercentCorrect*100);
        
        % This code makes plots that demonstrates how percent correct
        % depends on the criterion to which the likelihood ratio is
        % compared.  It should be best at 0, but isn't if you don't use the
        % ANALYTIC_LIKELY option above.
        DEBUG_CRITERIA = false;
        if (DEBUG_CRITERIA)
            llCriteria = linspace(min(logLikelyRatio),max(logLikelyRatio),1000);
            for jj = 1:length(llCriteria)
                for ii = 1:nTeObservations
                    if (logLikelyRatio(ii) > llCriteria(jj))
                        tePredictAsFunctionOfCriterion(ii,jj) = 1;
                    else
                        tePredictAsFunctionOfCriterion(ii,jj) = 0;
                    end
                end
                
                nCorrectAsFunctionOfCriterion(jj) = length(find(tePredictAsFunctionOfCriterion(:,jj) == teClasses));
                usePercentCorrectAsFunctionOfCriterion(jj) = nCorrectAsFunctionOfCriterion(jj)/length(teClasses);
            end
            
            figure(2); clf;
            subplot(1,2,1);
            hist(logLikelyRatio,40);
            subplot(1,2,2); hold on
            plot(llCriteria,usePercentCorrectAsFunctionOfCriterion,'r','LineWidth',4);
            plot([0 0],[0 1],'g','LineWidth',3);
            ylim([0 1]);
            drawnow;
        end
        
        % Optionally do Geisler analytic d-prime calculation, from Geisler 1984.
        % 
        % This assumes that the distribution a transformation of the
        % likelihood under each class is approximately normal, so we plot
        % histograms and look.  It's pretty close, and the Geisler numbers
        % are quite close to what we get empirically above.
        DO_GEISLER_VERSION = false;
        if (DO_GEISLER_VERSION)
            alphaMean = noStimData.noiseFreeIsomerizations(:)';
            betaMean = stimData.noiseFreeIsomerizations(:)';
            
            % Get analytic ideal observer dPrime and fraction correct
            [analyticFractionCorrect,analyticDPrime] = analyticPoissonIdealObserver(alphaMean,betaMean);
            
            % This next bit probably breaks for cases where there are 0
            % mean catches in the alpha data.
            for ii = 1:nTeObservations
                if (teClasses(ii) == 0)
                    alphaData = teData(ii,1:size(teData,2)/2);
                    betaData = teData(ii,size(teData,2)/2+1:end);
                else
                    betaData = teData(ii,1:size(teData,2)/2);
                    alphaData = teData(ii,size(teData,2)/2+1:end);
                end
                
                ZGivenAlpha(ii)= sum(alphaData .* log(betaMean./alphaMean));
                ZGivenBeta(ii)= sum(betaData .* log(betaMean./alphaMean));
            end
            
            meanZGivenBeta = mean(ZGivenBeta);
            meanZGivenAlpha = mean(ZGivenAlpha);
            varZGivenBeta = var(ZGivenBeta);
            varZGivenAlpha = var(ZGivenAlpha);
            dprimeEmpirical = (meanZGivenBeta-meanZGivenAlpha) / sqrt(0.5*(varZGivenBeta+varZGivenAlpha));
            
            figure(1); clf;
            subplot(1,3,1); hold on;
            [n,x] = hist(ZGivenAlpha,40);
            bar(x,n);
            pred = length(ZGivenAlpha)*normpdf(x,meanZGivenAlpha,sqrt(varZGivenAlpha))*(x(2)-x(1));
            plot(x,pred,'g','LineWidth',4);
            subplot(1,3,2); hold on;
            [n,x] = hist(ZGivenBeta,40);
            bar(x,n);
            pred = length(ZGivenBeta)*normpdf(x,meanZGivenBeta,sqrt(varZGivenBeta))*(x(2)-x(1));
            plot(x,pred,'g','LineWidth',4);
            
            % We can get the empirical ROC curve from the Z values and ask what
            % TAFC percent correct corresponds to those.
            criteria = linspace(min([ZGivenAlpha,ZGivenBeta]),max([ZGivenAlpha,ZGivenBeta]),1000);
            for iii = 1:length(criteria)
                empiricalHitRate(iii) = length(find(ZGivenBeta > criteria(iii))) / length(ZGivenBeta);
                empiricalFARate(iii) = length(find(ZGivenAlpha > criteria(iii))) / length(ZGivenBeta);
            end
            empiricalFractionCorrect = -trapz([1 empiricalFARate 0],[1 empiricalHitRate 0]);
            subplot(1,3,3);
            plot([1 empiricalFARate 0],[1 empiricalHitRate 0],'r','LineWidth',4);
            xlabel('False Alarm Rate');
            ylabel('Hit Rate');
            title('Empirical ROC Curve');
            axis('square'); axis([0 1 0 1]);
            drawnow;
            
            % Report percent correct
            fprintf('\td'': %2.2f (%2.2f), d'' percent correct %2.2f%%, percent correct from Z histo ROC %2.2f%%\n',...
                analyticDPrime,dprimeEmpirical,analyticFractionCorrect*100,empiricalFractionCorrect*100);
        end
       
        % Final newline
        fprintf('\n');
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
        useStdErr = std(percentCorrect)/sqrt(length(percentCorrect));

        % Report percent correct
        fprintf('\tPercent correct: %2.2f%%\n\n', usePercentCorrect*100);
        
        % No figure with this method
        h = [];  
        
    otherwise
        error('Unknown classification method');
end

end