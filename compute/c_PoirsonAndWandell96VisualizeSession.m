function c_PoirsonAndWandell96VisualizeSession()
% Conduct batch runs using the c_PoirsonAndWandel executive script
%
    % How many instances to generate
    nTrainingSamples = 512;
    
    % Select a subset, such as:
    spatialFrequency = 10;
    meanLuminance = 200;
    emPathType = 'frozen';
    freezeNoise = false;
    
    c_PoirsonAndWandell96Replicate(...
            'spatialFrequency', spatialFrequency, 'meanLuminance', meanLuminance, ...
            'nTrainingSamples', nTrainingSamples, 'emPathType', emPathType, ...
            'freezeNoise', freezeNoise, ...
            'computeResponses', false, 'visualizeResponses', true, ...
            'displayTrialBlockPartitionDiagnostics', false,  ...
            'findPerformance', false); 
end