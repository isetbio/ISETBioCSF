function [usePercentCorrect,useStdErr] = ClassifyForOneDirection(ii,data,theStimData,classes,nTrials,testContrasts,signalSource,nIntervals,PCAComponents,kFold)
% [usePercentCorrect,useStdErr] = ClassifyForOneDirection(ii,data,theStimData,classes,nTrials,testContrasts,signalSource,nIntervals,PCAComponents,kFold)
%
% Called within a parfor loop in t_colorGaborDetectFindPerformance.  We had
% to put this into a function to avoid some parfor errors.
%
% 7/14/16  dhb, xd  Pulled this out.  Hoping it will work.

responseSize = numel(theStimData{1}.responseInstanceArray(1).theMosaicPhotoCurrents(:));
for jj = 1:numel(testContrasts)
    fprintf('\nInserting (%d,%d) stimulus data from %d trials into design matrix ...', ii, jj, nTrials);
    if (nIntervals == 1)
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
    elseif (nIntervals == 2)
        for iTrial = 1:nTrials/2
            if (strcmp(signalSource,'photocurrents'))
                data(iTrial,responseSize+1:end) = theStimData{jj}.responseInstanceArray(iTrial).theMosaicPhotoCurrents(:);
                data(nTrials/2+iTrial,1:responseSize) = theStimData{jj}.responseInstanceArray(nTrials/2+iTrial).theMosaicPhotoCurrents(:);
            else
                data(iTrial,responseSize+1:end) = theStimData{jj}.responseInstanceArray(iTrial).theMosaicIsomerizations(:);
                data(nTrials/2+iTrial,1:responseSize) = theStimData{jj}.responseInstanceArray(nTrials/2+iTrial).theMosaicIsomerizations(:);
            end
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