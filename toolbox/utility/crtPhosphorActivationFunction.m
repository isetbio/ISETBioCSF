function phosphorFunction = crtPhosphorActivationFunction(refreshRate, samplesPerRefreshCycle) 
% phosphorFunction = crtPhosphorActivationFunction(refreshRate, samplesPerRefreshCycle) 
%
% Create a phoshor activarion function with sharp rise, shower decline
%
%  7/7/16  npc Wrote it.
%

    alpha = 1.9; t_50 = 0.02/1000; n = 2;
    phosphorFunction.timeInSeconds = linspace(0,1.0/refreshRate, samplesPerRefreshCycle);
    phosphorFunction.activation = (phosphorFunction.timeInSeconds.^n)./(phosphorFunction.timeInSeconds.^(alpha*n) + t_50^(alpha*n));
    phosphorFunction.activation = phosphorFunction.activation - phosphorFunction.activation(end);
    phosphorFunction.activation(phosphorFunction.activation<0) = 0;
    phosphorFunction.activation = phosphorFunction.activation / max(phosphorFunction.activation);
%     figure(3);
%     plot(phosphorFunction.timeInSeconds, phosphorFunction.activation, 'ko-');
%     pause
end

