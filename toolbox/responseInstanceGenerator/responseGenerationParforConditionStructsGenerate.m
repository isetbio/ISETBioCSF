function parforConditionStructs = responseGenerationParforConditionStructsGenerate(testConeContrasts,testContrasts)
% parforConditionStructs = responseGenerationParforConditionStructsGenerate(testConeContrasts,testContrasts)
%
% Create crossed list of contrast directions and contrasts to run and store
% in a struct array.

nParforConditions = size(testConeContrasts,2)*numel(testContrasts);
parforConditionStructs = cell(nParforConditions,1);
conditionIndex = 1;
for ii = 1:size(testConeContrasts,2) 
    for jj = 1:numel(testContrasts)
        thisConditionStruct.ii = ii;
        thisConditionStruct.jj = jj;
        thisConditionStruct.testConeContrasts = testConeContrasts(:,ii);
        thisConditionStruct.contrast = testContrasts(jj);
        parforConditionStructs{conditionIndex} = thisConditionStruct;
        conditionIndex = conditionIndex + 1;
    end
end