function testImagePipeline

    theRetina = makeRetina();
    
     imageSize = [128, 128, 3];
     input = ones(imageSize) * 1;
     [theRetina, ~, ~, ~, coneVec] = retinaCompute(theRetina, input);
%     
    
    %% Linearity test
    % Run N images with only one pixel set to a value of 1
    nInput = 10;
    coneResp = zeros(length(coneVec), nInput);
    oiMeanIlluminance2DAll = [];
    for idx = 1 : nInput
        input = zeros(imageSize);
        input(60+idx, 64,1) = 1;

        [theRetina,~, theOI, ~, coneVec] = retinaCompute(theRetina, input);
        oiMeanIlluminance2D = oiGet(theOI, 'illuminance');
        oiMeanIlluminance(idx) = mean(oiMeanIlluminance2D(:));
        if (isempty(oiMeanIlluminance2DAll))
            oiMeanIlluminance2DAll = oiMeanIlluminance2D;
        else
            oiMeanIlluminance2DAll = oiMeanIlluminance2DAll+oiMeanIlluminance2D;
        end
        % Save cone response for each image
        coneResp(:, idx) = coneVec;
    end

    % Run 1 image with N pixel set to the value of 1
    input = zeros(imageSize);
    input(60+(1:nInput), 64,1) = 1;
    [theRetina,~, theOI, ~, coneVec] = retinaCompute(theRetina, input);
    oiMeanIlluminance2D = oiGet(theOI, 'illuminance');
    oiMeanIlluminanceComposite = mean(oiMeanIlluminance2D(:));
    
    sum(oiMeanIlluminance) - oiMeanIlluminanceComposite
    figure(9999)
    imagesc(oiMeanIlluminance2D-oiMeanIlluminance2DAll)
    
    %% The sum of the cone response should be equal to the cone response to the sum of images
    % Scatter plot
    figure(); hold on;
    scatter(coneVec, sum(coneResp, 2));

    refPoint = [0, 20];
    plot(refPoint, refPoint);
    axis square;
    xlim(refPoint);
    ylim(refPoint);

end

function [theRetina, excitation, theOI, linearizedImage, allCone, L, M, S] = ...
    retinaCompute(theRetina, image)

    meanLuminanceCdPerM2 = [];
    [realizedStimulusScene, ~, linearizedImage] = sceneFromFile(image, 'rgb', ...
                meanLuminanceCdPerM2, theRetina.Display);
            

    % set the angular scene width
    realizedStimulusScene = sceneSet(realizedStimulusScene, 'fov', theRetina.FovealDegree);
            
    sceneMeanLuminance = sceneGet(realizedStimulusScene, 'luminance');
    sceneMeanLuminance = mean(sceneMeanLuminance(:));
    
    % optics
    theOI = oiCompute(theRetina.PSF, realizedStimulusScene);
    
    oiMeanIlluminance = oiGet(theOI, 'illuminance');
    oiMeanIlluminance = mean(oiMeanIlluminance(:));
    
    theRetina.LastOI = theOI;

    % cone excitations
    nTrialsNum = 2;
    emPath = zeros(nTrialsNum, 1, 2);

    % compute mosaic excitation responses
    % without the eye movement path
    theRetina.LastResponse = theRetina.Mosaic.compute(theOI, 'emPath', emPath);

    sizeExci   = size(theRetina.LastResponse);
    excitation = reshape(theRetina.LastResponse(1, :, :), [sizeExci(2), sizeExci(3)]);

    % LMS cone excitation
    L = getConetypeResponse(theRetina, theRetina.L_Cone_Idx);
    M = getConetypeResponse(theRetina, theRetina.M_Cone_Idx);
    S = getConetypeResponse(theRetina, theRetina.S_Cone_Idx);

    allCone = [L; M; S];
            
    figure()
    imagesc(sceneGet(realizedStimulusScene, 'RGB'));
    title(sprintf('mean lum: %g\n mean illum: %f\n', sceneMeanLuminance, oiMeanIlluminance));
    [sceneMeanLuminance oiMeanIlluminance]
    
    axis 'image'
    drawnow

end


function response = getConetypeResponse(this, type)
    sizeExci = size(this.LastResponse);
    excitation = reshape(this.LastResponse(1, :, :), [sizeExci(2), sizeExci(3)]);
    response   = excitation(this.Mosaic.pattern == type);
end
        
function this = makeRetina()
    this.FovealDegree = 1.0;
    eccBasedConeDensity = true;
    eccBasedConeQuantal = true;
    generateMosaic = ~true;
    if (generateMosaic)
        theMosaic = coneMosaicHex(5, ...                               % hex lattice sampling factor
                    'fovDegs',  this.FovealDegree , ...                           % match mosaic width to stimulus size
                    'eccBasedConeDensity', eccBasedConeDensity, ...            % cone density varies with eccentricity
                    'eccBasedConeQuantalEfficiency', eccBasedConeQuantal, ...  % cone quantal efficiency varies with eccentricity
                    'integrationTime', 0.2, ...                                % 0.1s integration time
                    'maxGridAdjustmentIterations', 50);  
             
        save('theMosaic.mat', 'theMosaic', '-v7.3');
        fprintf('Mosaic saved');
    else
        load('theMosaic.mat', 'theMosaic');
    end
    
    % Poisson noise model, mean response
    theMosaic.noiseFlag = 'none';
    this.Mosaic = theMosaic;
            
    % Display
    this.Display = displayCreate('LCD-Apple');
    this.Display.dist = 0.57;
            
    wavefrontSpatialSamples = 201;  % 501 results in 8.3262e-06
                                    % 901 results in 8.2817e-06
    pupilDiamMm = 3;
    umPerDegree = 300;
    
    % Point spread function of human eye
    this.PSF = oiCreate('wvf human', pupilDiamMm );
    % default results in illuminance difference of 3.3612e-05
                                    % 101 results in 7.0714e-06
                                    % 201 results in 8.1095e-06
    
    %this.PSF = oiWithCustomOptics('WvfHuman', wavefrontSpatialSamples, pupilDiamMm, umPerDegree);
       
    
    this.Mosaic.visualizeGrid(...
                'backgroundColor', [1 1 1], ...
                'ticksInVisualDegs', true);
            set(gca, 'XTick', [], 'YTick', []);
            set(gca,'YDir','reverse');
            
    this.L_Cone_Idx = 2;
    this.M_Cone_Idx = 3;
    this.S_Cone_Idx = 4;
end

