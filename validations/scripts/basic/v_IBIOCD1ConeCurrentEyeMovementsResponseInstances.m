function varargout = v_IBIOCD1ConeCurrentEyeMovementsResponseInstances(varargin)
% varargout = v_IBIOCD1ConeCuurentEyeMovementsResponseInstances(varargin)
%
% Works by running t_coneCurrentEyeMovementsResponseInstances with various arguments and comparing
% results with those stored.
%
% The 1 in the filename is to make sure that's gets run in the right order
% on a run through all validation scripts.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Hello
    UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_IBIOCD1ConeCurentEyeMovementsResponseInstances *****');
        %% Basic validation - no eye movements
    [validationData1, extraData1] = t_coneCurrentEyeMovementsResponseInstances(...
        'employStandardHostComputerResources', true, ...
        'generatePlots', runTimeParams.generatePlots, ...
        'visualizedResponseNormalization', 'submosaicBasedZscore', ...
        'freezeNoise', true);
    UnitTest.validationData('validationData1',validationData1);
    UnitTest.extraData('extraData1',extraData1);
    
    %% Spot version - no eye movements
    [validationData2, extraData2] = t_coneCurrentEyeMovementsResponseInstancesSpot(...
        'employStandardHostComputerResources', true, ...
        'generatePlots', runTimeParams.generatePlots,  ...
        'visualizedResponseNormalization', 'LMSabsoluteResponseBased', ...
        'freezeNoise', true);
    UnitTest.validationData('validationData2',validationData2);
    UnitTest.extraData('extraData2',extraData2);

    %% Photocurrent validation with frozen emPaths - rect mosaic
    rParams = responseParamsGenerate;
    rParams.spatialParams.gaussianFWHMDegs = 0.35;
    rParams.spatialParams.cyclesPerDegree = 8;
    rParams.spatialParams.fieldOfViewDegs = 0.6;
    rParams.temporalParams.secondsToInclude = 0.24;
    rParams.temporalParams.emPathType = 'frozen';
    rParams.mosaicParams.conePacking = 'rect';
    rParams.mosaicParams.fieldOfViewDegs = rParams.spatialParams.fieldOfViewDegs;
    rParams.mosaicParams.isomerizationNoise = 'frozen';
    rParams.mosaicParams.osNoise = 'frozen';

    % Select a small number of contitions
    testDirectionParams = instanceParamsGenerate();
    testDirectionParams.nAngles = 1;
    testDirectionParams.startAngle = 45;
    testDirectionParams.nContrastsPerDirection = 1;
    testDirectionParams.lowContrast = 0.9;
    testDirectionParams.highContrast = 0.9;
    testDirectionParams.trialsNum = 5;
    
    [validationData3, extraData3] = t_coneCurrentEyeMovementsResponseInstances(...
        'employStandardHostComputerResources', true, ...
        'rParams', rParams, ...
        'testDirectionParams', testDirectionParams, ...
        'generatePlots', runTimeParams.generatePlots, ...
        'visualizationFormat', 'montage', ...
        'visualizedResponseNormalization', 'submosaicBasedZscore', ...
        'freezeNoise', true);
    UnitTest.validationData('validationData3',validationData3);
    UnitTest.extraData('extraData3',extraData3);

    %% Photocurrent validation with frozen emPaths - hex mosaic
    rParams.mosaicParams.conePacking = 'hex';
    testDirectionParams.trialsNum = 2;

    [validationData4, extraData4] = t_coneCurrentEyeMovementsResponseInstances(...
        'employStandardHostComputerResources', true, ...
        'rParams', rParams, ...
        'testDirectionParams', testDirectionParams, ...
        'generatePlots', runTimeParams.generatePlots, ...
        'visualizationFormat', 'montage', ...
        'visualizedResponseNormalization', 'submosaicBasedZscore', ...
        'freezeNoise', true);
    UnitTest.validationData('validationData4',validationData4);
    UnitTest.extraData('extraData4',extraData4);
end




