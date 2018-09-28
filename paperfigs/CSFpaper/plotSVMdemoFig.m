function plotSVMdemoFig

    N = 100;
    rng default  % For reproducibility
    mu = [1 3];
    sigma = [4 1.0; 5.5 3];
    sigma = (sigma.' + sigma) / 2;
    
    features1 = mvnrnd(mu,sigma,N);
    class1Indices = 1:N;
    
    mu = [-3 3];
    sigma = [4 .1; 0.3 8];
    sigma = (sigma.' + sigma) / 2;
    features2 = mvnrnd(mu,sigma,N);
    class2Indices = N+(1:N);
    
    features = cat(1, features1, features2);
    
    features = bsxfun(@minus, features, mean(features,1));
    mean(features,1)
    
    classes(class1Indices) = 1;
    classes(class2Indices) = -1;
    svmStruct = svmtrain(features,classes,'ShowPlot',false, 'autoscale', false);

    
    sv = svmStruct.SupportVectors;
    alphaHat = svmStruct.Alpha;
    bias = svmStruct.Bias;
    kfun = svmStruct.KernelFunction;
    kfunargs = svmStruct.KernelFunctionArgs;
    N = 200;
    x = linspace(min(features(:,1)), max(features(:,1)), N);
    y = linspace(min(features(:,2)), max(features(:,2)), N);
	
    XYrange = min([max(x(:)) max(-x(:)) max(y(:)) max(-y(:))]) * [-1 1];
    %XYrange = 1.05*max([max(abs(x(:))) max(abs(y(:)))]) * [-1 1];
    
    [X,Y] = meshgrid(x, y);
    Xnew = [X(:) Y(:)];

    decisionBoundary = (feval(kfun,sv,Xnew,kfunargs{:})'*alphaHat(:)) + bias;
    z = max(abs(decisionBoundary(:)));

    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 600 600]);
    subplot('Position', [0.03 0.03 0.94 0.94]);
    imagesc(x,y,reshape(decisionBoundary,[N N]));
    hold on;
    [C,h] = contour(x,y,reshape(decisionBoundary,[N N]), [0 0]);
    h.LineColor = [0 0 0];
    h.LineWidth = 2.0;
    hold on
    plot(features(class1Indices,1), features(class1Indices,2), 'ko', 'MarkerSize', 18, 'MarkerFaceColor', [1 0.5 0.5]);
    plot(features(class2Indices,1), features(class2Indices,2), 'ko', 'MarkerSize', 18, 'MarkerFaceColor', [0.5 0.5 1]);
    %plot(sv(:,1), sv(:,2), 'ko', 'MarkerSize', 18);
    axis 'square'
    axis 'xy'
    box off;
    set(gca, 'XTick', [], 'YTick', []);
    set(gca, 'CLim', z*[-1 1], 'XLim', XYrange, 'YLim', XYrange);
    colormap(brewermap(1024, 'RdBu'));
    
    localDir = strrep(isetRootPath, 'toolboxes/isetbio/isettools', 'projects/IBIOColorDetect/paperfigs/CSFpaper/exports');
    NicePlot.exportFigToPDF(sprintf('%s/SVMdemo.pdf', localDir), hFig, 300);
end

