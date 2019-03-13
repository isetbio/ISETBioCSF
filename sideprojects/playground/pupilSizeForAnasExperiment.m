function pupilSizeForAnasExperiment
 
    % fraction of stimulus area to background (entire display) area
    fraction = 0.288;
    scene.xDegs = 20;
    scene.yDegs = 16.7;
    scene.luminance =  17.25;
    background.xDegs = scene.xDegs/sqrt(fraction);
    background.yDegs = scene.yDegs/sqrt(fraction);
    background.luminance = 0.364;
    
    exp1Params.age = 21; 
    exp1Params.eyeNum = 1; 
    exp1Params.area = background.xDegs*background.yDegs;
    exp1Params.adaptingFieldLuminance = computeMeanLuminance(scene, background);
    exp1Params.modelName = 'wy';  % Watson-Yelott
    pupilSizeForExp1 = humanPupilSize(exp1Params.adaptingFieldLuminance,exp1Params.modelName, exp1Params) 

    fraction = 0.3438;
    scene.xDegs = 18.56;
    scene.yDegs = 17.27;
    scene.luminance = 16.51;
    background.xDegs = scene.xDegs/sqrt(fraction);
    background.yDegs = scene.yDegs/sqrt(fraction);
    background.luminance = 0.43;
    % fraction of stimulus area to background (entire display) area
    
    exp2Params.age = 20; 
    exp2Params.eyeNum = 2; 
    exp2Params.area = background.xDegs*background.yDegs;
    exp2Params.adaptingFieldLuminance = computeMeanLuminance(scene, background);
    exp2Params.modelName = 'wy';  % Watson-Yelott
    pupilSizeForExp2Neutral = humanPupilSize(exp2Params.adaptingFieldLuminance,exp2Params.modelName, exp2Params) 

    scene.luminance = 15.29;
    exp2Params.adaptingFieldLuminance = computeMeanLuminance(scene, background);
    pupilSizeForExp2Red = humanPupilSize(exp2Params.adaptingFieldLuminance,exp2Params.modelName, exp2Params) 

    scene.luminance = 24.34;
    exp2Params.adaptingFieldLuminance = computeMeanLuminance(scene, background);
    pupilSizeForExp2Yellow = humanPupilSize(exp2Params.adaptingFieldLuminance,exp2Params.modelName, exp2Params) 

    
end

function meanLuminance = computeMeanLuminance(scene, background)
    l = zeros(round(background.xDegs*10), round(background.yDegs*10)) + background.luminance;
    l(1:round(scene.xDegs*10), 1:round(scene.yDegs*10)) = scene.luminance;
    meanLuminance = mean(l(:));
    figure()
    imagesc(l)
    axis 'image'
    set(gca, 'CLim', [0 30]);
    colormap(gray);
end
