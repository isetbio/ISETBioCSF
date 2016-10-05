function validationData = t_colorDetectFindPerformance(varargin)
% validationData = t_colorDetectFindPerformance(varargin)
%
% Classify data generated by
%   t_coneCurrentEyeMovementsResponseInstances.
% That tutorial generates multiple noisy instances of responses for color
% Gabors and saves them out.  Here we read the output and use an
% SVM to build a computational observer that gives us percent correct, and
% then do this for multiple contrasts.  The output of this tutorial can
% then be used by
%   t_colorGaborDetectThresholdsOnLMPlane
% to find and plot thresholds.
%
% See also
%   t_coneCurrentEyeMovementsResponseInstances,
%   t_colorGaborDetectIllustrateClassificationBoundary,
%   t_plotGabotDetectThresholdsOnLMPlane.
%
% Key/value pairs
%   'rParams' - Value the is the rParams structure to use
%   'testDirectionParams' - Value is the testDirectionParams structure to use
%   'setRngSeed' - true/false (default true).  Set the rng seed to a
%        value so output is reproducible.
%   'compute' - true/false (default true).  Do the computations.
%   'plotPsychometric' - true/false (default false).  Produce
%       psychometric function output graphs.
%   'plotSvmBoundary' - true/false (default false).  Plot classification boundary
%   'plotPCAAxis1' - First PCA component to plot (default 1)
%   'plotPCAAxis2' - Second PCA component to plot (default 2)
%   'delete' - true/false (default true).  Delete the response instance
%        files.  Useful for cleaning up utput when we are done with
%        it.  If this is true, output files are deleted at the end.
%        Does not currently delete figures or parent directories, but maybe
%        it should.

%% Parse input
p = inputParser;
p.addParameter('rParams',[],@isemptyorstruct);
p.addParameter('testDirectionParams',[],@isemptyorstruct);
p.addParameter('thresholdParams',[],@isemptyorstruct);
p.addParameter('setRng',true,@islogical);
p.addParameter('compute',true,@islogical);
p.addParameter('plotPsychometric',false,@islogical);
p.addParameter('plotSvmBoundary',false,@islogical);
p.addParameter('plotPCAAxis1',1,@isnumeric)
p.addParameter('plotPCAAxis2',2,@isnumeric)
p.addParameter('delete',false',@islogical);
p.parse(varargin{:});
rParams = p.Results.rParams;
testDirectionParams = p.Results.testDirectionParams;
thresholdParams = p.Results.thresholdParams;

%% Clear
if (nargin == 0)
    ieInit; close all;
end

%% Fix random number generator so we can validate output exactly
if (p.Results.setRng)
    rng(1);
end

%% Initialize validation data to empty
%
% It only gets set if we compute.
validationData = [];

%% Get the parameters we need
%
% t_colorGaborResponseGenerationParams returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
if (isempty(rParams))
    rParams = responseParamsGenerate;
    
    % Override some defult parameters
    %
    % Set duration equal to sampling interval to do just one frame.
    rParams.temporalParams.simulationTimeStepSecs = 200/1000;
    rParams.temporalParams.stimulusDurationInSeconds = rParams.temporalParams.simulationTimeStepSecs;
    rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.simulationTimeStepSecs;
    rParams.temporalParams.secondsToInclude = rParams.temporalParams.simulationTimeStepSecs;
    rParams.temporalParams.eyesDoNotMove = true;
    
    rParams.mosaicParams.timeStepInSeconds = rParams.temporalParams.simulationTimeStepSecs;
    rParams.mosaicParams.integrationTimeInSeconds = rParams.mosaicParams.timeStepInSeconds;
    rParams.mosaicParams.isomerizationNoise = true;
    rParams.mosaicParams.osNoise = true;
    rParams.mosaicParams.osModel = 'Linear';
end

%% Parameters that define the LM instances we'll generate here
%
% Make these numbers small (trialNum = 2, deltaAngle = 180,
% nContrastsPerDirection = 2) to run through a test quickly.
if (isempty(testDirectionParams))
    testDirectionParams = instanceParamsGenerate;
end

%% Parameters related to how we find thresholds from responses
if (isempty(thresholdParams))
    thresholdParams = thresholdParamsGenerate;
end

%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
readProgram = 't_coneCurrentEyeMovementsResponseInstances';
writeProgram = mfilename;

%% Compute if desired
if (p.Results.compute)
    
    % Read data for the no stimulus condition
    fprintf('Reading no stimulus data ... ');
    colorModulationParamsTemp = rParams.colorModulationParams;
    colorModulationParamsTemp.coneContrasts = [0 0 0]';
    colorModulationParamsTemp.contrast = 0;
    paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
    noStimData = rwObject.read('responseInstances',paramsList,readProgram);
    ancillaryData = rwObject.read('ancillaryData',paramsList,readProgram);
    
    % Get out some data we'll want
    nTrials = numel(noStimData.responseInstanceArray);
    testConeContrasts = ancillaryData.testConeContrasts;
    testContrasts = ancillaryData.testContrasts;
    
    % If everything is working right, these check parameter structures will
    % match what we used to specify the file we read in.
    %
    % SHOULD ACTUALLY CHECK FOR EQUALITY HERE.  Should be able to use
    % RecursivelyCompareStructs to do so.
    rParamsCheck = ancillaryData.rParams;
    LMPlaneInstanceParamsCheck = ancillaryData.LMPlaneInstanceParams;
    fprintf('done\n');
    
    % Do SVM for each test contrast and color direction.
    %
    % The work is done inside routine classifyForOneDirectionAndContrast.  We needed to
    % encapsulate it there to make parfor happy.
    %
    % If you don't have a computer configured to work with parfor, you may need
    % to change the parfor here to a plain for loop.
    tic
    parforConditionStructs = responseGenerationParforConditionStructsGenerate(testConeContrasts,testContrasts);
    nParforConditions = length(parforConditionStructs);
    usePercentCorrect = zeros(size(testConeContrasts,2),1);
    useStdErr = zeros(size(testConeContrasts,2),1);
   for kk = 1:nParforConditions
        thisConditionStruct = parforConditionStructs{kk};
        colorModulationParamsTemp = rParams.colorModulationParams;
        colorModulationParamsTemp.coneContrasts = thisConditionStruct.testConeContrasts;
        colorModulationParamsTemp.contrast = thisConditionStruct.contrast;
        paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp};
        stimData = rwObject.read('responseInstances',paramsList,readProgram);
        if (numel(stimData.responseInstanceArray) ~= nTrials)
            error('Inconsistent number of trials');
        end
        
        % Get performance for this condition.  Optional parameters control
        % whether or not the routine returns a handle to a plot that
        % illustrates the classifier.
        [usePercentCorrect(kk),useStdErr(kk),h] = ...
            classifyForOneDirectionAndContrast(noStimData,stimData,thresholdParams, ...
            'plotSvmBoundary',p.Results.plotSvmBoundary,'plotPCAAxis1',p.Results.plotPCAAxis1,'plotPCAAxis2',p.Results.plotPCAAxis2);
        
        % Save classifier plot if we made one and then close the figure.
        if (p.Results.plotSvmBoundary)
            paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, colorModulationParamsTemp, thresholdParams};
            rwObject.write(sprintf('svmBoundary_PCA%d_PCA%d',plotSvmPCAAxis1,plotSvmPCAAxis2), ...
                h,paramsList,writeProgram,'Type','figure');
            close(h);
        end
    end
    fprintf('Classification took %2.2f minutes\n', toc/60);
    clearvars('theData','useData','classificationData','classes');
    
    %% Take the returned vector form of the performance data and put it back into the
    % matrix form we expect below and elsewhere.
    %
    % See function responseGenerationParforConditionStructsGenerate for how we
    % pack the conditions into the order that this unpacks.
    for kk = 1:nParforConditions
        thisConditionStruct = parforConditionStructs{kk};
        performanceData.percentCorrect(thisConditionStruct.ii,thisConditionStruct.jj) = usePercentCorrect(kk);
        performanceData.stdErr(thisConditionStruct.ii,thisConditionStruct.jj) = useStdErr(kk);
    end
    
    %% Tuck away other information that we want to store
    performanceData.testConeContrasts = testConeContrasts;
    performanceData.testContrasts = testContrasts;
    performanceData.rParams = rParams;
    performanceData.LMPlaneInstanceParams = testDirectionParams;
    performanceData.thresholdParams = thresholdParams;
    clearvars('usePercentCorrect','useStdErr');
    
    %% Save classification performance data and a copy of this script
    fprintf('Writing performance data ... ');
    paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, thresholdParams};
    rwObject.write('performanceData',performanceData,paramsList,writeProgram);
    fprintf('done\n');
    
    %% Validation data
    if (nargin > 0)
        validationData = [];
    end
end

%% Plot performances obtained in each color direction as raw psychometric functions
if (p.Results.plotPsychometric) 
    fprintf('Reading performance data ... ');
    paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, thresholdParams};
    performanceData = rwObject.read('performanceData',paramsList,writeProgram);
    fprintf('done\n');
    
    for ii = 1:size(testConeContrasts,2)
        hFig = figure(1); clf;
        errorbar(testContrasts, squeeze(performanceData.percentCorrect(ii,:)), squeeze(performanceData.stdErr(ii, :)), ...
            'ro-', 'LineWidth', rParams.plotParams.lineWidth, 'MarkerSize', rParams.plotParams.markerSize, 'MarkerFaceColor', [1.0 0.5 0.50]);
        axis 'square'
        set(gca, 'YLim', [0 1.0],'XLim', [testContrasts(1) testContrasts(end)], 'FontSize', rParams.plotParams.axisFontSize);
        xlabel('contrast', 'FontSize' ,rParams.plotParams.labelFontSize, 'FontWeight', 'bold');
        ylabel('percent correct', 'FontSize' ,rParams.plotParams.labelFontSize, 'FontWeight', 'bold');
        box off; grid on
        title(sprintf('LMS = [%2.2f %2.2f %2.2f]', testConeContrasts(1,ii), testConeContrasts(2,ii), testConeContrasts(3,ii)), ...
            'FontSize',rParams.plotParams.titleFontSize);
        rwObject.write(sprintf('performanceData_%d',ii),hFig,paramsList,writeProgram,'Type','figure');
    end
end

%% Delete output data if desired
%
% Doesn't delete figures.  
if (p.Results.delete)
    paramsList = {rParams.spatialParams, rParams.temporalParams, rParams.oiParams, rParams.mosaicParams, rParams.backgroundParams, testDirectionParams, thresholdParams};
    rwObject.delete('performanceData',paramsList,writeProgram);
end

