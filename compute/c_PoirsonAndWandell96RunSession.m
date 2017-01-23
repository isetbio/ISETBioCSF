function c_PoirsonAndWandell96RunSession()
% Conduct batch runs using the c_PoirsonAndWandel executive script
%
    % Parameters varied
    nTrainingSamples = 512;
    spatialFrequency = 2;  %. 2,10
    meanLuminance = 200;  % 20,200
        
    % Actions to perform
    computeResponses = true;
    visualizeResponses = false;
    findPerformances = false;
    emPathType = 'random'; %select from 'random', 'frozen', 'frozen0', 'none'
    
    tBegin = clock;
    
    % Go!
    if (computeResponses)
        emPathType = 'frozen0';
        spatialFrequency = 2;  meanLuminance = 20;
        computeTheResponses(spatialFrequency, meanLuminance, nTrainingSamples, emPathType);
        spatialFrequency = 2;  meanLuminance = 200;
        computeTheResponses(spatialFrequency, meanLuminance, nTrainingSamples, emPathType);
        spatialFrequency = 10;  meanLuminance = 200;
        computeTheResponses(spatialFrequency, meanLuminance, nTrainingSamples, emPathType);
    end
    
    if (visualizeResponses)
        c_PoirsonAndWandell96Replicate(...
                'spatialFrequency', spatialFrequency, 'meanLuminance', meanLuminance, ...
                'nTrainingSamples', nTrainingSamples, 'emPathType', 'random', ...
                'computeResponses', false, 'visualizeResponses', true, ...
                'displayTrialBlockPartitionDiagnostics', false,  ...
                'findPerformance', false);
    end
    
    if (findPerformances)
        %findPerformancesForDifferentEvidenceIntegrationTimes(spatialFrequency, meanLuminance, nTrainingSamples);
        
        spatialFrequency = 2;  meanLuminance = 20;
        findThePerformances(spatialFrequency, meanLuminance, nTrainingSamples);
        
        spatialFrequency = 2;  meanLuminance = 200;
        findThePerformances(spatialFrequency, meanLuminance, nTrainingSamples);
        
        spatialFrequency = 10;  meanLuminance = 200;
        findThePerformances(spatialFrequency, meanLuminance, nTrainingSamples);
    end
    
    tEnd = clock;
    timeLapsed = etime(tEnd,tBegin);
    fprintf('BATCH JOB: Completed in %.2f hours. \n', timeLapsed/60/60);

    
end

function computeTheResponses(spatialFrequency, meanLuminance, nTrainingSamples, emPathType)   

    
    c_PoirsonAndWandell96Replicate(...
        'spatialFrequency', spatialFrequency, ...
        'meanLuminance', meanLuminance, ...
        'nTrainingSamples', nTrainingSamples, ...
        'computeResponses', true, ...
        'emPathType', emPathType, ...
        'visualizeResponses', false, ...
        'findPerformance', false);
    
end


function findPerformancesForDifferentEvidenceIntegrationTimes(spatialFrequency, meanLuminance, nTrainingSamples)

    emPathType = 'random';
    classifier = 'mlpt';
    performanceSignal = 'isomerizations';
    
    % Examine a range of integration times
    evidenceIntegrationTimes =   ([6 18 30 48 60 78 90 108 120 138 150 168 186 210]-1); % (5:10:250); %-([6 18 30 48 60 78 90 108 120 138 150 168 186 210]-1); % (5:10:250);
    for k = 1:numel(evidenceIntegrationTimes)
        evidenceIntegrationTime = evidenceIntegrationTimes(k);
        fprintf(2, 'Finding performance for ''%s'' EMpaths using an %s classifier operating on %2.1f milliseconds of the %s signals.\n', emPathType, classifier, evidenceIntegrationTime, performanceSignal);
        c_PoirsonAndWandell96Replicate(...
                'spatialFrequency', spatialFrequency, ...
                'meanLuminance', meanLuminance, ...
                'nTrainingSamples', nTrainingSamples, ...
                'computeResponses', false, ...
                'emPathType', emPathType, ...
                'visualizeResponses', false, ...
                'findPerformance', true, ...
                'performanceSignal', performanceSignal, ...
                'performanceClassifier', classifier, ...
                'performanceEvidenceIntegrationTime', evidenceIntegrationTime ....
                );
    end % k
    
    % And the the full time course
    c_PoirsonAndWandell96Replicate(...
                'spatialFrequency', spatialFrequency, ...
                'meanLuminance', meanLuminance, ...
                'nTrainingSamples', nTrainingSamples, ...
                'computeResponses', false, ...
                'emPathType', emPathType, ...
                'visualizeResponses', false, ...
                'findPerformance', true, ...
                'performanceSignal', performanceSignal, ...
                'performanceClassifier', classifier, ...
                'performanceEvidenceIntegrationTime', [] ....
                );
            
            
end

function findThePerformances(spatialFrequency, meanLuminance, nTrainingSamples)  
    for tableColumn = 3:3  % 1:4
        switch tableColumn
            case 1 
                % First column: mlpt on isomerizations 
                % for all 3 path types
                emPathTypes = {'none', 'random'};
                classifier = 'mlpt';
                performanceSignal = 'isomerizations';
            case 2
                % Second column: svm on isomerizations 
                % for all 3 path types
                emPathTypes = {'none', 'random'};
                classifier = 'svm';
                performanceSignal = 'isomerizations';
            case 3
                % Third column: svm on photocurrents 
                % for the 2 non-static path types
                emPathTypes = {'random'};
                classifier = 'svm';
                performanceSignal = 'photocurrents';
            case 4
                % Fourth column: svm on V1 filter bank (operating on photocurrents) 
                % for the 2 non-static path types
                emPathTypes = {'random'};
                classifier = 'svmV1FilterBank';
                performanceSignal = 'photocurrents';
            otherwise
                error('No params for table column: %d', tableColumn);
        end % switch

        for emPathTypeIndex = 1:numel(emPathTypes)
            emPathType = emPathTypes{emPathTypeIndex};
            fprintf(2, 'Finding performance for ''%s'' EMpaths using an %s classifier on the %s signals.\n', emPathType, classifier, performanceSignal);
            c_PoirsonAndWandell96Replicate(...
                'spatialFrequency', spatialFrequency, ...
                'meanLuminance', meanLuminance, ...
                'nTrainingSamples', nTrainingSamples, ...
                'computeResponses', false, ...
                'emPathType', emPathType, ...
                'visualizeResponses', false, ...
                'findPerformance', true, ...
                'performanceSignal', performanceSignal, ...
                'performanceClassifier', classifier ...
                );
        end %  for emPathTypeIndex
    end % tableColumn
 end
