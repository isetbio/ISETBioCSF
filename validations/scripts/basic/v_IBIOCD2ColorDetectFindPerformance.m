function varargout = v_IBIOCD2ColorDetectFindPerformance(varargin)
% varargout = v_IBIOCDC2olorDetectFindPerformance(varargin)
%
% Works by running t_colorDetectFindPerformance with various arguments and comparing
% results with those stored.
%
% The 2 in the filename is to make sure that's gets run in the right order
% on a run through all validation scripts.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Hello
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOCD2ColorDetectFindPerformance *****');
    
        %% Get run params
    rParams = responseParamsGenerate;
    
    % Override some defult parameters
    rParams.spatialParams.fieldOfViewDegs = 1.0;
    rParams.spatialParams.cyclesPerDegree = 8;
    rParams.spatialParams.gaussianFWHMDegs = 0.375;
    % Set duration equal to sampling interval to do just one frame.
    rParams.temporalParams.stimulusDurationInSeconds = 200/1000;
    rParams.temporalParams.stimulusSamplingIntervalInSeconds = rParams.temporalParams.stimulusDurationInSeconds;
    rParams.temporalParams.secondsToInclude = rParams.temporalParams.stimulusDurationInSeconds;
    % Mosaic params
    rParams.mosaicParams.integrationTimeInSeconds = rParams.temporalParams.stimulusDurationInSeconds;
    rParams.mosaicParams.isomerizationNoise = 'frozen';         % Type coneMosaic.validNoiseFlags to get valid values
    rParams.mosaicParams.osNoise = 'frozen';                    % Type outerSegment.validNoiseFlags to get valid values
    rParams.mosaicParams.osModel = 'Linear';
    rParams.mosaicParams.fieldOfViewDegs = 1.0;
    
    testDirectionParams = instanceParamsGenerate;
    testDirectionParams.trialsNum = 64;
    
    %% Basic validation
    [validationData1,extraData1] = t_colorDetectFindPerformance(...
        'rParams', rParams, ...
        'employStandardHostComputerResources', true,...
        'generatePlots', runTimeParams.generatePlots,...
        'plotPsychometric', true,...
        'freezeNoise',true);
    UnitTest.validationData('validationData1',validationData1, ...
        'UsingTheFollowingVariableTolerancePairs', ...
         'validationData1.percentCorrect', 0.11, ...
         'validationData1.stdErr', 0.06 ...
        );
    UnitTest.extraData('extraData1',extraData1);
end