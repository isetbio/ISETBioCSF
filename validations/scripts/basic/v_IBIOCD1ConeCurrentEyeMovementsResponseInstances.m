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
        'generatePlots', runTimeParams.generatePlots, ...
        'freezeNoise', true);
    UnitTest.validationData('validationData1',validationData1);
    UnitTest.extraData('extraData1',extraData1);
    
    %% Spot version - no eye movements
    [validationData2, extraData2] = t_coneCurrentEyeMovementsResponseInstancesSpot(...
        'generatePlots', runTimeParams.generatePlots,  ...
        'freezeNoise', true);
    UnitTest.validationData('validationData2',validationData2);
    UnitTest.extraData('extraData2',extraData2);
    
    %% Photocurrent validation with frozen emPaths
    %
    % Use a mosaic that covers the central 1/3 of the stimulus to accelerate computation
    rParams = responseParamsGenerate;
    rParams.temporalParams.secondsToInclude = 0.5;
    rParams.mosaicParams.fieldOfViewDegs = 0.3*rParams.spatialParams.fieldOfViewDegs;
    rParams.mosaicParams.isomerizationNoise = 'frozen';
    rParams.mosaicParams.osNoise = 'frozen';
    
    % Select a small number of contitions
    testDirectionParams = instanceParamsGenerate();
    testDirectionParams.nAngles = 1;
    testDirectionParams.nContrastsPerDirection = 1;
    testDirectionParams.lowContrast = 0.9;
    testDirectionParams.highContrast = 0.9;
    testDirectionParams.trialsNum = 20;
    
    [validationData3, extraData3] = t_coneCurrentEyeMovementsResponseInstances(...
        'rParams', rParams, ...
        'testDirectionParams', testDirectionParams, ...
        'emPathType', 'frozen', ...
        'generatePlots', runTimeParams.generatePlots, ...
        'freezeNoise', true);
    UnitTest.validationData('validationData3',validationData3);
    UnitTest.extraData('extraData3',extraData3);


end




