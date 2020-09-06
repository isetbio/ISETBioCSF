function plotSVMdemoFig

    N = 100;
    rng default  % For reproducibility
    mu = [1 3];
    sigma = [4 1.0; 5.5 3];
    sigma = (sigma.' + sigma) / 2;
    
    % Generate data
    features1 = mvnrnd(mu,sigma,N);
    class1Indices = 1:N;
    
    % Add some covariance
    mu = [-3 3];
    sigma = [4 .1; 0.3 8];
    sigma = (sigma.' + sigma) / 2;
    features2 = mvnrnd(mu,sigma,N);
    class2Indices = N+(1:N);
    
    % Setup features vector
    features = cat(1, features1, features2);
    features = bsxfun(@minus, features, mean(features,1));
    
    % Setup classes vector
    classes(class1Indices) = 1;
    classes(class2Indices) = -1;

    % Fit an SVM
    svm = fitcsvm(features,classes, ...
                'KernelFunction', 'linear', ...
                'Standardize', true ...
                );

    %decisionBoundaryFunction = @(x) -(x*svm.Beta(1) + svm.Bias)/svm.Beta(2);
    %sv = svm.SupportVectors;

    % Generate high-resolution xy grid  for plotting the boundary
    N = 200;
    x = linspace(min(features(:,1)), max(features(:,1)), N);
    y = linspace(min(features(:,2)), max(features(:,2)), N);
    [X,Y] = meshgrid(x, y);

    % Run the SVM to get assigned classes for all the xy grid points
    pred = [X(:),Y(:)];
    [p, score] = predict(svm,pred);
    
    % Get our decision boundary (score)
    decisionBoundary = score(:,1);
   
    % Render the plot using plotlab
    renderPlot(x, y, decisionBoundary, features, class1Indices, class2Indices);
end

function renderPlot(x,y, decisionBoundary, features, class1Indices, class2Indices)
    
    % Setup figure
    figWidthInches = 10; figHeightInches = 10;
    plotlabOBJ = setupPlotLab(0, figWidthInches, figHeightInches);
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1]);
    
    % Generate axes
    ax = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', 1, ...
        'colsNum', 1, ...
        'leftMargin', 0.01, ...
        'widthMargin', 0.03, ...
        'heightMargin', 0.04, ...
        'bottomMargin', 0.03, ...
        'rightMargin', 0.01, ...
        'topMargin', 0.03);
    
    ax = ax{1,1};
    N = length(x);
    XYrange = min([max(x(:)) max(-x(:)) max(y(:)) max(-y(:))]) * [-1 1];
    z = max(abs(decisionBoundary(:)));

    % The decision boundary as a density plot
    imagesc(ax,x,y,reshape(decisionBoundary,[N N]));
    hold(ax, 'on');
    
    % The decision boundary as a line
    [C,h] = contour(ax,x,y,reshape(decisionBoundary,[N N]), [0 0]);
    h.LineColor = [0 0 0];
    h.LineWidth = 2.0;
    
    % The data points
    scatter(ax,features(class1Indices,1), features(class1Indices,2));
    scatter(ax,features(class2Indices,1), features(class2Indices,2));
    
    % Finalize figure
    axis(ax, 'square'); axis(ax, 'xy'); box(ax,'on');
    set(ax, 'XTick', [], 'YTick', []);
    set(ax, 'CLim', z*[-1 1], 'XLim', XYrange, 'YLim', XYrange);
    colormap(ax,brewermap(1024, 'RdYlGn'));
    
    % Export to PDF
    localDir = strrep(isetRootPath, 'toolboxes/isetbio/isettools', 'projects/ISETBioCSF/paperfigs/CSFpaper/exports');
    NicePlot.exportFigToPDF(sprintf('%s/SVMdemo.pdf', localDir), hFig, 300);
end


function plotlabOBJ = setupPlotLab(mode, figWidthInches, figHeightInches)
    if (mode == 0)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [0.8 0.1 0.1; 0.1 0.5 0.1], ...
                'axesBox', 'off', ...
                'axesTickDir', 'both', ...
                'renderer', 'painters', ...
                'lineWidth', 6, ...
                'lineMarkerSize', 24, ...
                'contourLineWidth', 6, ...
                'axesFontSize', 30, ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'figureWidthInches', figWidthInches, ...
                'figureHeightInches', figHeightInches);
    else
        pause(2.0);
        plotlab.resetAllDefaults();
    end
end 

