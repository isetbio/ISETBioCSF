function c_PoirsonAndWandellReplicate
% c_PoirsonAndWandellReplicate(varargin)
%
% Compute color detection thresholds to replicate the Poirson & Wandell 1996

% Start with default parameters
rParams = responseParamsGenerate;
LMplaneDirectionParams = instanceParamsGenerate;
thresholdParams = thresholdParamsGenerate;

% Adapt spatial params to match P&W 1996
displayFlattenedStruct = true;
rParams = adaptSpatialParamsBasedOnConstantCycleCondition(rParams, displayFlattenedStruct);

% Adapt LMplane params to match P&W 1996
LMplaneDirectionParams = adaptLMPLaneDirectionParamsBasedOnFig3A(LMplaneDirectionParams, displayFlattenedStruct);       

% Start computing
theOI = colorDetectOpticalImageConstruct(rParams.oiParams)

end


% ---- HELPER FUNCTIONS ----
function LMplaneDirectionParams = adaptLMPLaneDirectionParamsBasedOnFig3A(LMplaneDirectionParams, displayFlattenedStruct)

    % Adapt to the Figure 3A params of P&W 1996 params
    LMplaneDirectionParams = setExistingFieldsInStruct(LMplaneDirectionParams, ...
        {              'startAngle',  45; ...
                       'deltaAngle',  90; ...
                          'nAngles',  1; ...
            'nContrastsPerDirection', 20; ...        % Number of contrasts to run in each color direction
                      'lowContrast',  0.0001; ...
                     'highContrast',  0.1;
                    'contrastScale', 'log' ...      % Choose between 'linear' and 'log'
        });
    
    if (displayFlattenedStruct)
        LMplaneDirectionParamsFlattened = UnitTest.displayNicelyFormattedStruct(LMplaneDirectionParams, 'LMplaneDirectionParams', '', 65) 
    end
end

function rParams = adaptSpatialParamsBasedOnConstantCycleCondition(rParams, displayFlattenedStruct)
    % Adapt to the Figure 3 params of P&W 1996 params
    % - constant cycle condition:  Gaussian spatial window decreases as spatial frequency increases
    % - SF = 2 cpd, 
    % - width of the Gaussian window at half height = 1.9 deg
    % viewing distance: 0.75 meters
    rParams.spatialParams = setExistingFieldsInStruct(rParams.spatialParams, ...
        {    'cyclesPerDegree', 2.0; ...
            'gaussianFWHMDegs', 1.9; ...
             'fieldOfViewDegs', 5.0; ...        % In P&W 1996, in the constan cycle condition, this was 10 deg (Section 2.2, p 517)
             'viewingDistance', 0.75; ...
                         'row', 256; ...        % stim height in pixels
                         'col', 256 ...         % stim width in pixels
        });

    
    % In the constant cycle condition, the background was xyY= 0.38, 0.39, 536.2 cd/m2
    % Also they say they placed a uniform field to the screen to increase the contrast resolutionso (page 517, Experimental Aparratus section)
    rParams.backgroundParams = setExistingFieldsInStruct(rParams.backgroundParams, ...
    	{   'backgroundxyY', [0.38 0.39 536.2]; ... 
        	  'monitorFile', 'CRT-MODEL'; ...
               'leakageLum', 0.5 ...                    % I (Nicolas) suspect this is too low
    	});
    
    % In the constant cycle condition the display (Barco CDCT6351) was refreshed at 87 Hz
    CRTrefreshRateInHz = 87;
    rParams.temporalParams = setExistingFieldsInStruct(rParams.temporalParams, ...
        {                   'frameRate', CRTrefreshRateInHz; ...
               'simulationTimeStepSecs', 1.0/CRTrefreshRateInHz; ...    % scompute responses once per CRT refresh cycle; this will need to bearound 0.1 milliseconds for biophys-based outersegment
                   'windowTauInSeconds', 0.165; ...                     % 165 msec
    'stimulusSamplingIntervalInSeconds', 1.0/CRTrefreshRateInHz; ...        % sample the stimulus once per CRT refresh cycle
            'stimulusDurationInSeconds', 5*0.165;  ...      % 5 * window Tau
                     'secondsToInclude', 1.0; ...           % analyze the cental +/- 0.5 seconds of the response
               'secondsToIncludeOffset', 0.0; ...           % within the stimulus peak
                        'eyesDoNotMove', true, ...          % let's assume no eye movements for now
        });
    
    % Update computed temporal properties
    [sampleTimes,gaussianTemporalWindow, rasterModulation] = gaussianTemporalWindowCreate(rParams.temporalParams);
    rParams.temporalParams = setExistingFieldsInStruct(rParams.temporalParams, ...
        {           'sampleTimes', sampleTimes; ...
         'gaussianTemporalWindow', gaussianTemporalWindow; ...
                   'nSampleTimes', length(sampleTimes) ...
        });
    
    
    % Adapt optical image params
    rParams.oiParams = setExistingFieldsInStruct(rParams.oiParams, ...
        {   'fieldOfViewDegs', rParams.spatialParams.fieldOfViewDegs*1.5 ...  % make it's FOV  50% larger than the stimulus
        });
    
    
    % Adapt mosaic params 
    rParams.mosaicParams = setExistingFieldsInStruct(rParams.mosaicParams, ...
        {         'fieldOfViewDegs', rParams.spatialParams.fieldOfViewDegs; ...           % match stimulus size
                'timeStepInSeconds', rParams.temporalParams.simulationTimeStepSecs; ...
         'integrationTimeInSeconds', rParams.temporalParams.simulationTimeStepSecs; ...
               'isomerizationNoise', true; ...
                          'osNoise', true; ...
                          'osModel', 'Linear'...
        });
    
    
    if (displayFlattenedStruct)
        rParamsFlattened = UnitTest.displayNicelyFormattedStruct(rParams, 'rParams', '', 65) 
    end
    
end



