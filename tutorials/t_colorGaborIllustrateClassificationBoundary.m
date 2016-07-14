%% t_colorGaborIllustrateClassificationBoundary
%
% Plots the data and svm boundary for a specified condition projected onto
% the first 2 principal components.
%
% 7/13/16  xd  wrote it based on code in t_colorGaborDetectFindThresholds

%% Initialize
ieInit; clear; close all;

% Add project toolbox to Matlab path
AddToMatlabPathDynamically(fullfile(fileparts(which(mfilename)),'../toolbox'));

%% Define parameters of analysis
%
% Condition directory that has the response instances
conditionDir = 'cpd2_sfv1.00_fw0.350_tau0.165_dur0.33_em0_use50_off35_b0_l1_LMS1.00_0.00_0.00_mfv1.00';

% Signal source: select between 'photocurrents' and 'isomerizations'
signalSource = 'photocurrents';

% Number of SVM cross validations to use
kFold = 5;

% PCA components.  Set to zero for no PCA
PCAComponents = 200;

% Specify the conditions we would like to plot. 
ColorDirection = 6;
contrastLevel = 5;

%% Get data saved by t_colorGaborConeCurrentEyeMovementsResponseInstances
%
% The file ending in _0 gives us the no stimulus data
dataDir = colorGaborDetectOutputDir(conditionDir,'output');
responseFile = 'responseInstances_0';
responsesFullFile = fullfile(dataDir, sprintf('%s.mat',responseFile));
classificationPerformanceFile = fullfile(dataDir, sprintf('ClassificationPerformance_%s_kFold%0.0f_pca%0.0f.mat',signalSource,kFold,PCAComponents));
fprintf('\nLoading data from %s ...', responsesFullFile);
theBlankData = load(responsesFullFile);
fprintf('done\n');
nTrials = numel(theBlankData.theNoStimData.responseInstanceArray);
testContrasts = theBlankData.testContrasts;

%% Put zero contrast response instances into data that we will pass to the SVM
responseSize = numel(theBlankData.theNoStimData.responseInstanceArray(1).theMosaicPhotoCurrents(:));
fprintf('\nInserting null stimulus data from %d trials into design matrix ...\n', nTrials);
for iTrial = 1:nTrials
    if (iTrial == 1)
        data = zeros(2*nTrials, responseSize);
        classes = zeros(2*nTrials, 1);
    end
    if (strcmp(signalSource,'photocurrents'))
        data(iTrial,:) = theBlankData.theNoStimData.responseInstanceArray(iTrial).theMosaicPhotoCurrents(:);
    else
        data(iTrial,:) = theBlankData.theNoStimData.responseInstanceArray(iTrial).theMosaicIsomerizations(:);
    end
    
    % Set up classes variable
    classes(iTrial,1) = 0;
    classes(nTrials+iTrial,1) = 1;
end
fprintf('done\n');

%% Do SVM for each test contrast and color direction.
%
% The work is done inside routine ClassifyForOneDirection.  We needed to
% encapsulate it there to make parfor happy.
%
% If you don't have a computer configured to work with parfor, you may need
% to change the parfor just to plain for.
tic
thisResponseFile = sprintf('responseInstances_%d',ColorDirection);
thisResponseFullFile = fullfile(dataDir, sprintf('%s.mat',thisResponseFile));
theStimData = load(thisResponseFullFile);
theStimData = theStimData.theStimData;

fprintf('\nInserting (%d,%d) stimulus data from %d trials into design matrix ...', ColorDirection, contrastLevel, nTrials);
for iTrial = 1:nTrials
    % Put data into the right form for SVM.
    % This loop overwrites the stimlus data each time through, a
    % little risky, coding wise, but will work unless someone
    % modifies the data generation tutorial to produce a different
    % number of noisy instances for different test directions or
    % contrasts.
    if (strcmp(signalSource,'photocurrents'))
        data(nTrials+iTrial,:) = theStimData{contrastLevel}.responseInstanceArray(iTrial).theMosaicPhotoCurrents(:);
    else
        data(nTrials+iTrial,:) = theStimData{contrastLevel}.responseInstanceArray(iTrial).theMosaicIsomerizations(:);
    end
end
fprintf(' done\n');

% Do PCA?
if (PCAComponents > 0)
    fprintf('\tDoing PCA ... ');
    theData = transformDataWithPCA(data,PCAComponents);
    fprintf('done\n');
else
    theData = data;
end

% Perform SVM classification for this stimulus vs the zero contrast stimulus
fprintf('\tRunning SVM for chromatic direction %d, contrast %2.2f ...',ColorDirection,testContrasts(contrastLevel));
[usePercentCorrect,useStdErr,svm] = classifyWithSVM(theData,classes,kFold);
fprintf(' correct: %2.2f%%\n', usePercentCorrect*100);
fprintf('SVM classification took %2.2f minutes\n', toc/60);

%% Plotting here

% Choose two principal components to serve as the axis of the figure
PC1 = 1;
PC2 = 100;

% Plot the data points along these two components
figure('Position',[0 0 750 750]);
hold on;
plot(theData(1:500,PC1),theData(1:500,PC2),'*','MarkerSize',7);
plot(theData(501:end,PC1),theData(501:end,PC2),'o','MarkerSize',7);

% Extract the boundary line. svm.Beta is the line orthogonal to the
% boundary. We take it's null space (and since we only have 2 dimensions,
% it is a single vector) which is the boundary. We scale it by an
% arbitrarily large number so the line comes out nicely.
orthLineToBoundary = svm.Beta;
orthLineToBoundary = orthLineToBoundary([PC1 PC2]);
boundary = null(orthLineToBoundary')*100;

% Plot the decision boundary
plot([-boundary(1) boundary(1)],[-boundary(2) boundary(2)],'--k','LineWidth',2);

% Reset the axis limits
ylim([min(theData(:,PC2)) max(theData(:,PC2))]);
xlim([min(theData(:,PC1)) max(theData(:,PC1))]);

% Make plot look nice + various labels/titles/legend
set(gca,'FontName','Helvetica','FontSize',18,'LineWidth',2);
legend({'Null Stimulus','Stimulus'},'Location','Northeast','FontSize',18);

axis square;
xlabel(sprintf('Principal Component %d',PC1),'FontSize',22);
ylabel(sprintf('Principal Component %d',PC2),'FontSize',22);
title('SVM Decision Boundary on Two Principal Components','FontSize',24);
