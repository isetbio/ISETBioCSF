function varargout = v_BanksEtAlReplicate(varargin)
% varargout = v_anksEtAlReplicate(varargin)
%
% Works by running t_coneIsomerizationsMovie with various arguments and comparing
% results with those stored.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Checks
%
% 1) - Quantal catch.  12/27/16.
%
% Geisler 1984 gives the formula he uses to estimate mean cone quantal
% catch.  For 100 ms and a 2 mm pupil, this yields 3571 quanta per cone at
% 340 cd/m2, for his cone aperture of 0.6 min (0.28 min2).  See  routine
% IsomerizationsFromLuminanceGeisler.
% 
% If we go look by hand at the cone catches for these parameters produced
% for the uniform field by t_coneCurrentEyeMovementsResponseInstances for
% the regular hex mosaic, we find 2956 for the M cones and 4097 for the L
% cones.  At 2:1 L:M, this gives an average of 3716, which is pretty darn
% close.)
%
% 2) - Pupil size.  12/27/16
% 
% Doubling the pupil size increase retinal irradiance by a factor of 4,
% which should in turn double contrast sensitivity.  For 75% correct
% threshold, 340 cd/m2, 10 cpd and with optical blur, sensitivity goes from
% 600 for a 2 mm pupil to 1268 for a 4 mm pupil, which seems close enough
% for a numerical calculation like this one although why it consistently
% seems to be less than half bothers me a little. This is for the hex
% mosaic.
%
% 3) - Changing criterion percent correct. 12/27/16.
%
% It's not clear whether Banks et al. used 75% correct as their ideal
% observer criterion (as Geisler did in his 1984 paper) or switched to
% 70.1% to match the 2 down - 1 up reversal threshold from the
% psychophysics.  Changing the criterion changes the sensitivity in the
% expected direction from the pupil size case just above: 2mm pupil -> 750;
% 4mm pupil -> 1610.  I am not quite sure why increasing the pupil size
% helps a tad more than a factor of 2 both here and above.
%
% 4) - Taking out blur.  12/27/16.  For the 2 mm pupil regular hex mosaic
% case, taking out the optical blur for case 3 above takes sensitivity to
% 1006, from 750, which is the correct direction for 10 cpd.

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
        freezeNoise = true;
    else
        doSimulationWithBanksEtAlmosaicParams = false;
        computeResponses = true;
        findPerformance = true;
        fitPsychometric = true;
        doBlur = true;
        nTrainingSamples = 500;
        thresholdCriterionFraction = 0.75;
        freezeNoise = false;
    end
    
    %% Basic validation
    if (doSimulationWithBanksEtAlmosaicParams)
        % Run with the Banks mosaic
        [validationData1, extraData1] = c_BanksEtAlReplicate('useScratchTopLevelDirName',true, ...
            'computeResponses',computeResponses,'findPerformance',findPerformance,'fitPsychometric',fitPsychometric,...
            'nTrainingSamples',nTrainingSamples,'thresholdCriterionFraction',thresholdCriterionFraction,...
            'conePacking','hexReg','innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(3),'coneSpacingMicrons', 3.0, ... 
            'blur',doBlur,'cyclesPerDegree',10,'luminances',340,'pupilDiamMm',2,'generatePlots',runTimeParams.generatePlots,'freezeNoise',freezeNoise);
        UnitTest.validationData('validationData1',validationData1);
        UnitTest.extraData('extraData1',extraData1);

        [validationData2, extraData2] = c_BanksEtAlReplicate('useScratchTopLevelDirName',true, ...
            'computeResponses',computeResponses,'findPerformance',findPerformance,'fitPsychometric',fitPsychometric,...
            'nTrainingSamples',nTrainingSamples,'thresholdCriterionFraction',thresholdCriterionFraction,...
            'conePacking','hexReg','innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(3),'coneSpacingMicrons', 3.0, ...
            'blur',doBlur,'cyclesPerDegree',10,'luminances',340,'pupilDiamMm',4,'generatePlots',runTimeParams.generatePlots,'freezeNoise',freezeNoise);
        UnitTest.validationData('validationData2',validationData2);
        UnitTest.extraData('extraData2',extraData2);
    else
        % Run with rect mosaic and parameters we think match September 2016.
        %   Rectangular mosaic, 2uM pixel size, 1.4uM collecting area.
        %
        %   Pupil size was 3mm because of the bug where it did not
        %   actually get changed.
        %
        % These parameters produce 1836 mean isomerizations for the M cones
        % and 2542 for the L cones, for the uniform background, and these
        % numbers are the same in December 2016 as they were in September
        % 2016.  The mean isomerizations are now represented at single
        % rather than double precision, however.  That in turn seems to
        % make no difference when I change it back now.
        % 
        % For 5000 training samples, constrast sensitivity is 532 in
        % December 2016.  This is lower than it was in September 2016
        c_BanksEtAlReplicate('useScratchTopLevelDirName',true, ...
            'computeResponses',computeResponses,'findPerformance',findPerformance,'fitPsychometric',fitPsychometric,...
            'nTrainingSamples',nTrainingSamples,'thresholdCriterionFraction',thresholdCriterionFraction,...
            'conePacking','rect','innerSegmentSizeMicrons',1.4,'coneSpacingMicrons',2, ...   
            'blur',doBlur,'cyclesPerDegree',10,'luminances',340,'pupilDiamMm',3,'generatePlots',runTimeParams.generatePlots,'freezeNoise',freezeNoise);
    end
    
end

