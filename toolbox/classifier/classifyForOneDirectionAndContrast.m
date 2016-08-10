function [usePercentCorrect,useStdErr,h] = classifyForOneDirectionAndContrast(stimData,classificationData,classes,thresholdParams,varargin)
% [usePercentCorrect,useStdErr,h] = classifyForOneDirectionAndContrast(stimData,data,classes,thresholdParams,varargin)
%
% Do classification for one stimulus color direction and contrast.
%
% Optional key/value pairs
%   'Plot' - true/false (default false).  Plot classification boundary
%   'PlotAxis1' - First PCA component to plot (default 1)
%   'PlotAxis2' - Second PCA component to plot (default 2)
%
% 7/14/16  dhb, xd  Pulled this out.  Hoping it will work.

% Parse args
p = inputParser;
p.addRequired('stimData',@isstruct);
p.addRequired('classificationData',@isnumeric);
p.addRequired('classes',@isnumeric);
p.addRequired('thresholdParams',@isstruct);
p.addParameter('Plot',false,@islogical);
p.addParameter('PlotAxis1',1,@isnumeric)
p.addParameter('PlotAxis2',2,@isnumeric)
p.parse(stimData,classificationData,classes,thresholdParams,varargin{:});

% Insert stimulus responses into data for classifier
[classificationData,classes] = classificationDataStimDataInsert(classificationData,classes,stimData,thresholdParams);

% Do PCA?
if (thresholdParams.PCAComponents > 0)
    fprintf('\tDoing PCA ... ');
    theData = transformDataWithPCA(classificationData,thresholdParams.PCAComponents);
    fprintf('done\n');
else
    theData = classificationData;
end

% Perform SVM classification for this stimulus vs the zero contrast stimulus
fprintf('\tRunning SVM ...');
[usePercentCorrect, useStdErr, svm] = classifyWithSVM(theData,classes,thresholdParams.kFold);
fprintf(' correct: %2.2f%%\n', usePercentCorrect*100);

% Optional plot of classification boundary
if (p.Results.Plot)
    % Plot the data points along these two components
    PCAAxis1 = p.Results.PlotAxis1;
    PCAAxis2 = p.Results.PlotAxis2;
    
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

end