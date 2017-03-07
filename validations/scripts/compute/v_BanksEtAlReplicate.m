function varargout = v_BanksEtAlReplicate(varargin)
% varargout = v_BanksEtAlReplicate(varargin)
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
% 4mm pupil -> 1610.  Again, I am not quite sure why increasing the pupil size
% helps a tad more than a factor of 2 both here and above.
%
% 4) - Taking out blur.  12/27/16.  For the 2 mm pupil regular hex mosaic
% case, taking out the optical blur for case 3 above takes sensitivity to
% 1006, from 750, which is the correct direction for 10 cpd.
%
% The numbers above are before I switched the maximum likelihood classifier
% calculation from using poisspdf to a calculation that takes advantage of
% the analytic form of the Poisson pdf. Set ANALYTIC_LIKELY to false in
% classifyForOneDirectionAndContrast to get the old behavior.  It is now
% set to true.  I have also subsequently reduced the number of contrasts
% that the routine computes for -- go back to defaults to produce old
% behavior exactly.
%
% Also reduced wavelength sampling and number of training samples for speed
% in the validation Original values are [380 4 780] and 100, for comparison
% with numbers above.  The wavelength sampling actually matters in terms of
% the threshold you get.  Finally changed to 20 cpd from 10 cpd because
% this is a smaller number of cones and thus also makes things run faster.

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    %% Hello
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOCDIBanksEtAlReplicate *****');
    
    %% Freeze reandom number generator seed
    rng('default');
    
    %% Parameters
    VALIDATE = true;
    if (VALIDATE)
        doSimulationWithBanksEtAlmosaicParams = true;
        params.computeResponses = true;
        params.findPerformance = true;
        params.fitPsychometric = true;
        params.LMSRatio = [0.62 0.31 0.07];
        params.blur = true;
        params.apertureBlur = false;
        params.nTrainingSamples = 50;
        params.thresholdCriterionFraction = 0.75;
        params.freezeNoise = true;
        params.useTrialBlocks = false;
        params.wavelengths = [400 20 700];
    
    % This runs with old parameters.  Probably don't need it anymore, as I
    % added this case to the main validation file.  The notes may be of
    % interest someday if we ever need to back way out.
    else
        doSimulationWithBanksEtAlmosaicParams = false;
        params.computeResponses = true;
        params.findPerformance = true;
        params.fitPsychometric = true;
        params.LMSRatio = [0.62 0.31 0.07];
        params.blur = true;
        params.apertureBlur = false;
        params.nTrainingSamples = 500;
        params.thresholdCriterionFraction = 0.75;
        params.freezeNoise = false;
        params.useTrialBlocks = false;

    end
    
    %% Basic validation
    %
    % This branch checks the regular hex mosaic with 2 pupil sizes
    % and one rect mosaic, and doesn't check against validation data.
    if (doSimulationWithBanksEtAlmosaicParams)
        % Run with the Banks mosaic
        [validationData1, extraData1] = c_BanksEtAlReplicate(params, ...
            'useScratchTopLevelDirName',true, 'generatePlots',runTimeParams.generatePlots, ...
            'conePacking','hexReg','innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(3),'coneSpacingMicrons', 3.0, ...
            'cyclesPerDegree',20,'luminances',340,'pupilDiamMm',2,...
            'nContrastsPerDirection',4,'lowContrast',0.005,'highContrast',0.1,'contrastScale','log');
        UnitTest.validationData('validationData1',validationData1, ...
            'UsingTheFollowingVariableTolerancePairs', ...
            'validationData1.mlptThresholds.thresholdContrasts', 1e-3);
        UnitTest.extraData('extraData1',extraData1);

        [validationData2, extraData2] = c_BanksEtAlReplicate(params, ...
            'useScratchTopLevelDirName',true, 'generatePlots',runTimeParams.generatePlots, ...
            'conePacking','hexReg','innerSegmentSizeMicrons', sizeForSquareApertureFromDiameterForCircularAperture(3),'coneSpacingMicrons', 3.0, ...
            'cyclesPerDegree',20,'luminances',340,'pupilDiamMm',4,...
            'nContrastsPerDirection',4,'lowContrast',0.005,'highContrast',0.1,'contrastScale','log');
        validationValue = validationData2.mlptThresholds.thresholdContrasts;
        UnitTest.validationData('validationData2',validationValue, ...
            'UsingTheFollowingVariableTolerancePairs', ...
            'validationValue', 1e-6);
        UnitTest.extraData('extraData2',extraData2);
        
        [validationData3, extraData3] = c_BanksEtAlReplicate(params, ...
            'useScratchTopLevelDirName',true,'generatePlots',runTimeParams.generatePlots, ...
            'conePacking','rect','innerSegmentSizeMicrons',1.4,'coneSpacingMicrons',2, ...
            'cyclesPerDegree',20,'luminances',340,'pupilDiamMm',3, ...
            'nContrastsPerDirection',4,'lowContrast',0.005,'highContrast',0.1,'contrastScale','log');
        validationValue = validationData3.mlptThresholds.thresholdContrasts;
        UnitTest.validationData('validationData3',validationValue, ...
            'UsingTheFollowingVariableTolerancePairs', ...
            'validationValue', 1e-6);
        UnitTest.extraData('extraData3',extraData3);
        
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
        % 2016.
        % 
        % For 5000 training samples, constrast sensitivity is 532 in
        % December 2016.  This is lower than it was in September 2016,
        % where it was about 726.  This difference turns out to be because
        % we now compute in single rather than double precision, and the
        % log likelihood signal known exactly classifier actually does a
        % bit better with a double precision representation of the mean
        % isomerizations.  Along the way to finding this out, I verified
        % that:
        %
        % Although the noise free isomerizations are very close to same between
        % September and December, when we use the September code to
        % classify December response instances, we get the December
        % performance.  Thus the difference is in the 
        % responses.  The difference is subtle - nothing pops out from
        % histogramming the responses to the L cones for the uniform field
        % or for a grating, and the noise free isomerizations are equal to
        % about 1 part in 10^8.
        %
        % The noise is not frozen, and the difference over time is real
        %  September thresholds (%): 0.0624, 0.0630, 0.0671, 0.0651
        %  December thresholds (%):  0.0884, 0.0996, 0.0912
        %  December double prec (%): 0.0627
        % The variability is apparent in the individual points on the
        % psychometric function, as well as in the thresholds.
        %
        % Both September and December code drops to 50% correct at low
        % constasts, but this drop happens at lower contrasts (below
        % 0.0001) for the September code.  Performance for both times
        % points is 50.00% for a 0% contrast stimulus, which it should
        % because then there really is no difference between the two
        % classes being classified.
        %
        % The numbers above are before I switched the maximum likelihood
        % classifier calculation from using poisspdf to a calculation that takes
        % advantage of the analytic form of the Poisson pdf. Set ANALYTIC_LIKELY to
        % false in classifyForOneDirectionAndContrast to get the old behavior.  It
        % is now set to true.
        c_BanksEtAlReplicate(params, ...
            'useScratchTopLevelDirName',true,'generatePlots',runTimeParams.generatePlots, ...
            'conePacking','rect','innerSegmentSizeMicrons',1.4,'coneSpacingMicrons',2, ...   
            'cyclesPerDegree',10,'luminances',340,'pupilDiamMm',3);
    end
    
end

