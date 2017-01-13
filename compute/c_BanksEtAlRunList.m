% c_BanksEtAlRunList
%
% Run the Banks et al. computations with various set parameters.

% Clear and close
clear; close all;

% Common parameters
params.useScratchTopLevelDirName = false;
params.computeResponses = false;
params.findPerformance = false;
params.fitPsychometric = false;
params.nTrainingSamples = 2000;
params.conePacking = 'hexReg';
params.mosaicRotationDegs = 30;
params.coneSpacingMicrons = 3;
params.innerSegmentSizeMicrons = sizeForSquareApertureFromDiameterForCircularAperture(params.coneSpacingMicrons);
params.coneDarkNoiseRate = [0 0 0];
params.LMSRatio = [0.67 0.33 0];
params.cyclesPerDegree = [5 10 20 30 40 50 60];
params.luminances = [3.4 34 340];
params.pupilDiamMm = 2;
params.thresholdCriterionFraction = 0.701;
params.freezeNoise = true;
params.useTrialBlocks = true;
params.generatePlots = true;
params.plotCSF = true;
params.visualizeResponses = false;

% No blur, no dark noise, no aperture blur
params.blur = false;
params.apertureBlur = false;
params.coneDarkNoiseRate = [0 0 0];
c_BanksEtAlReplicate(params);

% No blur, no dark noise, with aperture blur
params.blur = false;
params.apertureBlur = true;
params.coneDarkNoiseRate = [0 0 0];
c_BanksEtAlReplicate(params);

% With blur, no dark noise, with aperture blur
params.blur = true;
params.apertureBlur = true;
params.coneDarkNoiseRate = [0 0 0];
c_BanksEtAlReplicate(params);

% % Plot Banks et al. data
% plotParams = plotParamsGenerate;
% [A,B,C,D,E] = LoadDigitizedBanksFigure2;
% 
% hFig = figure; clf; hold on
% fontBump = 4;
% markerBump = -4;
% lineBump = -1;
% set(gcf,'Position',[100 100 450 650]);
% set(gca,'FontSize', plotParams.axisFontSize+fontBump);
% plot(A(:,1),A(:,2),'r-','MarkerSize',plotParams.markerSize+markerBump,'LineWidth',plotParams.lineWidth+lineBump);
% plot(B(:,1),B(:,2),'r:','MarkerSize',plotParams.markerSize+markerBump,'LineWidth',plotParams.lineWidth+lineBump);
% plot(C(:,1),C(:,2),'b-','MarkerSize',plotParams.markerSize+markerBump,'LineWidth',plotParams.lineWidth+lineBump);
% plot(D(:,1),D(:,2),'b-','MarkerSize',plotParams.markerSize+markerBump,'MarkerFaceColor','r','LineWidth',plotParams.lineWidth+lineBump);
% plot(E(:,1),E(:,2),'b-','MarkerSize',plotParams.markerSize+markerBump,'MarkerFaceColor','r','LineWidth',plotParams.lineWidth+lineBump);
% xlabel('Log10 Spatial Frequency (cpd)', 'FontSize' ,plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
% ylabel('Log10 Contrast Sensitivity', 'FontSize' ,plotParams.labelFontSize+fontBump, 'FontWeight', 'bold');
% xlim([1 100]); ylim([10 10000]);
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% box off; grid on
% titleStr1 = 'Banks Et Al. CSFs';
% titleStr2 = '';
% title({titleStr1 ; titleStr2});


