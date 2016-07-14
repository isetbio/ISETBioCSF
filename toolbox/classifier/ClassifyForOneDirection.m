function [usePercentCorrect,useStdErr] = ClassifyForOneDirection(ii,data,theStimData,classes,nTrials,testContrasts,signalSource,PCAComponents,kFold)
% [usePercentCorrect,useStdErr] = ClassifyForOneDirection(ii,data,theStimData,classes,nTrials,testContrasts,signalSource,PCAComponents,kFold)
%
% Called within a parfor loop in t_colorGaborDetectFindThresholds.  We had
% to put this into a function to avoid some parfor errors.
%
% 7/14/16  dhb, xd  Pulled this out.  Hoping it will work.

for jj = 1:numel(testContrasts)
    fprintf('\nInserting (%d,%d) stimulus data from %d trials into design matrix ...', ii, jj, nTrials);
    for iTrial = 1:nTrials
        % Put data into the right form for SVM.
        % This loop overwrites the stimlus data each time through, a
        % little risky, coding wise, but will work unless someone
        % modifies the data generation tutorial to produce a different
        % number of noisy instances for different test directions or
        % contrasts.
        if (strcmp(signalSource,'photocurrents'))
            data(nTrials+iTrial,:) = theStimData{jj}.responseInstanceArray(iTrial).theMosaicPhotoCurrents(:);
        else
            data(nTrials+iTrial,:) = theStimData{jj}.responseInstanceArray(iTrial).theMosaicIsomerizations(:);
        end
    end
    fprintf(' done\n');
    
    % Do PCA?
    if (PCAComponents > 0)
        fprintf('\tDoing PCA ... ');
        theData = transformDataWithPCA(data,PCAComponents);
        fprintf('done\n');
    else
        theData = data;
    end
    
    % Perform SVM classification for this stimulus vs the zero contrast stimulus
    fprintf('\tRunning SVM for chromatic direction %d, contrast %2.2f ...', ii , testContrasts(jj));
    [usePercentCorrect(jj), useStdErr(jj)] = classifyWithSVM(theData,classes,kFold);
    fprintf(' correct: %2.2f%%\n', usePercentCorrect(jj)*100);
end

end