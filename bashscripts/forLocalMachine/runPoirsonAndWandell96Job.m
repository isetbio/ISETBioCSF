function runPoirsonAndWandell96Job()

    %tbUseProject('IBIOColorDetect');
    nTrials = 1024;
    performanceClassifierTrainingSamples = 768-64;  %1024-128-32
	%c_PoirsonAndWandell96RunSession(1, nTrials, performanceClassifierTrainingSamples);
    %c_PoirsonAndWandell96RunSession(2, nTrials, performanceClassifierTrainingSamples);
    c_PoirsonAndWandell96RunSession(3, nTrials, performanceClassifierTrainingSamples);
    c_PoirsonAndWandell96RunSession(4, nTrials, performanceClassifierTrainingSamples);
    c_PoirsonAndWandell96RunSession(5, nTrials, performanceClassifierTrainingSamples);
    c_PoirsonAndWandell96RunSession(6, nTrials, performanceClassifierTrainingSamples);

end

