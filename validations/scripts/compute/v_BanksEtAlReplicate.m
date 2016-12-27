function varargout = v_BanksEtAlReplicate(varargin)
% varargout = v_anksEtAlReplicate(varargin)
%
% Works by running t_coneIsomerizationsMovie with various arguments and comparing
% results with those stored.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    %% Hello
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOCDIBanksEtAlReplicate *****');
    
    %% Freeze reandom number generator seed
    rng('default');
    
    %% Parameters
    VALIDATE = false;
    if (VALIDATE)
        doSimulationWithBanksEtAlmosaicParams = true;
        computeResponses = true;
        findPerformance = true;
        fitPsychometric = true;
        doBlur = true;
        nTrainingSamples = 100;
        thresholdCriterionFraction = 0.75;
    else
        doSimulationWithBanksEtAlmosaicParams = false;
        computeResponses = true;
        findPerformance = true;
        fitPsychometric = true;
        doBlur = true;
        nTrainingSamples = 500;
        thresholdCriterionFraction = 0.75;
    end
    
    %% Basic validation
    if (doSimulationWithBanksEtAlmosaicParams)
        % Run with the Banks mosaic
        [validationData1, extraData1] = c_BanksEtAlReplicate('useScratchTopLevelDirName',true, ...
            'computeResponses',computeResponses,'findPerformance',findPerformance,'fitPsychometric',fitPsychometric,...
            'nTrainingSamples',nTrainingSamples,'thresholdCriterionFraction',thresholdCriterionFraction,...
            'conePacking','hexReg','innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(3),'coneSpacingMicrons', 3.0, ... 
            'blur',doBlur,'cyclesPerDegree',10,'luminances',340,'pupilDiamMm',2,'generatePlots',runTimeParams.generatePlots);
        UnitTest.validationData('validationData1',validationData1);
        UnitTest.extraData('extraData1',extraData1);

        [validationData2, extraData2] = c_BanksEtAlReplicate('useScratchTopLevelDirName',true, ...
            'computeResponses',computeResponses,'findPerformance',findPerformance,'fitPsychometric',fitPsychometric,...
            'nTrainingSamples',nTrainingSamples,'thresholdCriterionFraction',thresholdCriterionFraction,...
            'conePacking','hexReg','innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(3),'coneSpacingMicrons', 3.0, ...
            'blur',doBlur,'cyclesPerDegree',10,'luminances',340,'pupilDiamMm',4,'generatePlots',runTimeParams.generatePlots);
        UnitTest.validationData('validationData2',validationData2);
        UnitTest.extraData('extraData2',extraData2);
 
    else
        % Run with rect mosaic and parameters we think match summmer 2016.
        % 
        % Rectangular mosaic, 2uM pixel size, 1.4uM collecting area, we
        % think.
        %
        % And pupil size was 3mm because of the bug where it did not
        % actually get changed.
        c_BanksEtAlReplicate('useScratchTopLevelDirName',true, ...
            'computeResponses',computeResponses,'findPerformance',findPerformance,'fitPsychometric',fitPsychometric,...
            'nTrainingSamples',nTrainingSamples,'thresholdCriterionFraction',thresholdCriterionFraction,...
            'conePacking','rect','innerSegmentSizeMicrons',1.4,'coneSpacingMicrons',2, ...   
            'blur',doBlur,'cyclesPerDegree',10,'luminances',340,'pupilDiamMm',3,'generatePlots',runTimeParams.generatePlots);
    end
    
end

