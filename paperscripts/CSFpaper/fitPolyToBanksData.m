function fitPolyToBanksData

    load('BanksSubjects.mat', 'pjbSubjectSFs', 'pjbSubjectCSFs', 'msbSubjectSFs', 'msbSubjectCSFs');

    figure(1); clf;
    plot(pjbSubjectSFs, pjbSubjectCSFs, 'ro', 'MarkerSize', 12); hold on;
    plot(msbSubjectSFs, msbSubjectCSFs, 'bo', 'MarkerSize', 12);
    set(gca, 'XLim', [1 60], 'YLim', [1 200], 'XScale', 'log', 'YScale', 'log');
    

    [pjbSubjectSFsHiRes, pjbSubjectCSFsHiRes] = fitData(pjbSubjectSFs, pjbSubjectCSFs);
    plot(pjbSubjectSFsHiRes, pjbSubjectCSFsHiRes, 'r-', 'LineWidth', 1.5);
    
    [msbSubjectSFsHiRes, msbSubjectCSFsHiRes] = fitData(msbSubjectSFs, msbSubjectCSFs);
    plot(msbSubjectSFsHiRes, msbSubjectCSFsHiRes, 'b-', 'LineWidth', 1.5);
    
    meanSubjectCSF = 0.5*(msbSubjectCSFsHiRes+pjbSubjectCSFsHiRes);
    plot(msbSubjectSFsHiRes, meanSubjectCSF, 'k-', 'LineWidth', 1.5);
    targetSF = [2 4 8 16 32];
    for k = 1:numel(targetSF)
        idx = find(msbSubjectSFsHiRes == targetSF(k));
        if ~isnan(idx)
            plot(msbSubjectSFsHiRes(idx), meanSubjectCSF(idx), 'ks');
        end
    end
    
end

function [xFit, yFit] = fitData(xdata,ydata)

    xdata = xdata(:);
    ydata = ydata(:);
    
    F = @(p,xdata)p(1) + p(2)*exp(-p(3)*xdata.^p(6)) + p(4)*exp(-p(5)*xdata.^2);
    params0 = [1.3090   43.5452    0.0379   92.5870    0.0356    1.3392]; 
    
    [params,resnorm,~,exitflag,output] = lsqcurvefit(F,params0,xdata,ydata);
    
    params
    xFit = 2:0.1:30;
    yFit = F(params,xFit);
    
    
   % yFit = interp1(x(:),y(:), xFit, 'pchip');

end
