function [V1filterBank, hFig] = generateV1FilterBank(spatialParams, mosaicParams, topLevelDirParams, visualizeSpatialScheme, thresholdParams, paramsList)
    
    % Filter width
    filterWidthInDegrees = spatialParams.fieldOfViewDegs;
    if (isfield(thresholdParams,'poolingTemplateWidthInDegrees'))
        filterWidthInDegrees = thresholdParams.poolingTemplateWidthInDegrees;
        spatialParams.fieldOfViewDegs = filterWidthInDegrees;
    end
    
    % Load the mosaic
    coneParamsList = {topLevelDirParams, mosaicParams};
    theProgram = 't_coneCurrentEyeMovementsResponseInstances';
    rwObject = IBIOColorDetectReadWriteBasic;
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    
    % Get the cone locs in degrees
    if (strcmp(mosaicParams.conePacking,'rect'))
        coneLocsInMeters = theMosaic.coneLocs;
    else
        coneLocsInMeters = theMosaic.coneLocsHexGridAlignedWithSerializedConeMosaicResponse(); % OLD, incorrect way: coneLocsHexGrid;
    end
    coneLocsInDegs(:,1) = coneLocsInMeters(:,1) / theMosaic.width * theMosaic.fov(1);
    coneLocsInDegs(:,2) = coneLocsInMeters(:,2) / theMosaic.height * theMosaic.fov(2);
       
    if (strcmp(mosaicParams.conePacking, 'hex'))
        p = properties(theMosaic);
        if (any(ismember(p, 'eccBasedConeQuantalEfficiency'))) && (theMosaic.eccBasedConeQuantalEfficiency == false)
            % Compute density map around each cone, so we can normalize the
            % spatial pooling weights at each location based on the cone 
            % density at that location
            eccInMeters = sqrt(sum(coneLocsInMeters.^2, 2));
            ang = atan2(squeeze(coneLocsInMeters(:,2)), squeeze(coneLocsInMeters(:,1)))/pi*180;
            [~, ~, coneDensity] = coneSizeReadData('eccentricity',eccInMeters(:),'angle',ang(:));
        else
            % Cone efficiency corrections with eccentricity are on, so do
            % not normalize spatial pooling weights with respect to cone density
            coneDensity = 1;
        end
    else
        coneDensity = 1;
    end
    
    switch(spatialParams.spatialType)
        case 'Gabor'  
            % Make the stimulus spatial modulation
            spatialPattern = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, 1.0));
            spatialModulation = spatialPattern-1;
            spatialModulation = spatialModulation/max(abs(spatialModulation(:)));
            
        otherwise
            error('Currently generating V1 filter banks for Gabor stimuli only.');
    end  % switch
    
    xaxis = (0:(size(spatialModulation,2)-1))/size(spatialModulation,2) * abs(filterWidthInDegrees);
    xaxis = xaxis - mean(xaxis);
    yaxis = (0:(size(spatialModulation,1)-1))/size(spatialModulation,1) * abs(filterWidthInDegrees);
    yaxis = yaxis - mean(yaxis);
            
    % Generate the V1 filter bank        
    V1filterBank = makeV1FilterBank(spatialParams, filterWidthInDegrees, coneLocsInDegs, xaxis, yaxis, coneDensity, thresholdParams.spatialPoolingKernelParams);

    if (visualizeSpatialScheme)
        coneRadiusMicrons = 0.5*diameterForCircularApertureFromWidthForSquareAperture(theMosaic.pigment.width)*1e6;
        iTheta = (0:30:360) / 180 * pi;
        if (theMosaic.shouldCorrectAbsorptionsWithEccentricity())
            % Compute ecc-varying apertures
            dx = coneRadiusMicrons * 1e-6;
            apertureOutline.x = dx * cos(iTheta);
            apertureOutline.y = dx * sin(iTheta);
        
            [coneApertureOutline, ~, ~] = theMosaic.computeApertureSizes(...
                dx, [], apertureOutline, [], coneLocsInMeters(:,1), coneLocsInMeters(:,2));
            
            coneApertureOutline.x = coneApertureOutline.x *1e6/theMosaic.micronsPerDegree;
            coneApertureOutline.y = coneApertureOutline.y *1e6/theMosaic.micronsPerDegree;
        else
            coneRadiusDegs = coneRadiusMicrons/theMosaic.micronsPerDegree;
            coneApertureOutline.x = coneRadiusDegs * cos(iTheta);
            coneApertureOutline.y = coneRadiusDegs * sin(iTheta);
        end
        hFig = visualizeSpatialPoolingScheme(xaxis, yaxis, spatialModulation, ...
            thresholdParams.spatialPoolingKernelParams, V1filterBank, coneLocsInDegs, ...
            mosaicParams.fieldOfViewDegs, spatialParams.fieldOfViewDegs, coneApertureOutline);
    
        % Save figure
        theProgram = mfilename;
        rwObject = IBIOColorDetectReadWriteBasic;
        data = 0;
        fileName = sprintf('V1poolingKernel');
        paramsList{numel(paramsList)+1} = thresholdParams;
        rwObject.write(fileName, data, paramsList, theProgram, ...
           'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
       
    end % visualizeSpatialScheme      
end


function V1filterBank = makeV1FilterBank(spatialParams, filterWidthDegs, coneLocsDegs, xaxisDegs, yaxisDegs, coneDensity, spatialPoolingKernelParams)

    validV1KernelTypes = {'V1envelope', 'V1SinUnit', 'V1CosUnit', 'V1QuadraturePair'};
    if (~ismember(spatialPoolingKernelParams.type, validV1KernelTypes))
        validV1KernelTypes
        error('spatialPoolingKernelParams.type (''%s'') is not a valid V1 kernel type.\n', spatialPoolingKernelParams.type);
    end
    
    % filter width
    spatialParams.gaussianFWHMDegs = spatialPoolingKernelParams.shrinkageFactor * filterWidthDegs/2.0;
    
    % make the cos-phase filter
    spatialParams.ph = 0;
    cosPhaseFilter = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, 1.0));
    V1filterBank.cosPhasePoolingProfile = cosPhaseFilter-1;
    
    % make the sin-phase filter
    spatialParams.ph = pi/2;
    sinPhaseFilter = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, 1.0));
    V1filterBank.sinPhasePoolingProfile = sinPhaseFilter-1;

    
    % Compute energy envelope
    RFprofile = V1filterBank.cosPhasePoolingProfile.^2 + V1filterBank.sinPhasePoolingProfile.^2;
    V1filterBank.RFprofile = RFprofile / max(abs(RFprofile(:)));
    V1filterBank.xaxisDegs = xaxisDegs;
    V1filterBank.yaxisDegs = yaxisDegs;
    
    % Find pooling weights
    [X,Y] = meshgrid(xaxisDegs, yaxisDegs);
    rfCoordsDegs = [X(:) Y(:)];
    
    [~, idx] = pdist2(rfCoordsDegs, coneLocsDegs, 'euclidean', 'Smallest', 1);
    V1filterBank.cosPhasePoolingWeights = V1filterBank.cosPhasePoolingProfile(idx);
    V1filterBank.sinPhasePoolingWeights = V1filterBank.sinPhasePoolingProfile(idx);
    V1filterBank.envelopePoolingWeights = V1filterBank.RFprofile(idx);
    
    maxCos = max(abs(V1filterBank.cosPhasePoolingWeights(:)));
    maxSin = max(abs(V1filterBank.sinPhasePoolingWeights(:)));
    maxEnvelope = max(abs(V1filterBank.envelopePoolingWeights));
    
    % Adjust weights by the inverse of the coneDensity
    if (spatialPoolingKernelParams.adjustForConeDensity)
        V1filterBank.cosPhasePoolingWeights = V1filterBank.cosPhasePoolingWeights ./ coneDensity;
        V1filterBank.sinPhasePoolingWeights = V1filterBank.sinPhasePoolingWeights ./ coneDensity;
        V1filterBank.envelopePoolingWeights = V1filterBank.envelopePoolingWeights ./ coneDensity;
    end
    
    V1filterBank.cosPhasePoolingWeights = V1filterBank.cosPhasePoolingWeights / max(abs(V1filterBank.cosPhasePoolingWeights(:))) * maxCos;
    V1filterBank.sinPhasePoolingWeights = V1filterBank.sinPhasePoolingWeights / max(abs(V1filterBank.sinPhasePoolingWeights(:))) * maxSin;
    V1filterBank.envelopePoolingWeights = V1filterBank.envelopePoolingWeights / max(abs(V1filterBank.envelopePoolingWeights(:))) * maxEnvelope;
    
    % Normalize with respect to total energy over space
    netWeight = 1/sqrt(2.0) * sqrt(sum(V1filterBank.cosPhasePoolingWeights(:).^2) + sum(V1filterBank.sinPhasePoolingWeights(:).^2));
    V1filterBank.cosPhasePoolingWeights = V1filterBank.cosPhasePoolingWeights/netWeight;
    V1filterBank.sinPhasePoolingWeights = V1filterBank.sinPhasePoolingWeights/netWeight;
    
    if (any(V1filterBank.cosPhasePoolingWeights(:)<-1))
        fprintf('Cos pooling contains % values< -1\n', sum(V1filterBank.cosPhasePoolingWeights(:)<-1));
    end
    
    if (any(V1filterBank.cosPhasePoolingWeights(:)>1))
        fprintf('Cos pooling contains % values< >1\n', sum(V1filterBank.cosPhasePoolingWeights(:)>1));
    end
    
    if (any(V1filterBank.sinPhasePoolingWeights(:)<-1))
        fprintf('Sin pooling contains % values< -1\n', sum(V1filterBank.sinPhasePoolingWeights(:)<-1));
    end
    
    if (any(V1filterBank.sinPhasePoolingWeights(:)>1))
        fprintf('Sin pooling contains % values< >1\n', sum(V1filterBank.sinPhasePoolingWeights(:)>1));
    end
    
   
    if (strcmp(spatialPoolingKernelParams.type, 'V1SinUnit'))
        V1filterBank.cosPhasePoolingProfile = 0*V1filterBank.cosPhasePoolingProfile;
    end
    
    if (strcmp(spatialPoolingKernelParams.type, 'V1CosUnit'))
        V1filterBank.sinPhasePoolingProfile = 0*V1filterBank.sinPhasePoolingProfile;
    end
    
    if (any(isnan(V1filterBank.cosPhasePoolingWeights)))
        fprintf('Cos pooling contains %d NANs\n', sum(isnan(V1filterBank.cosPhasePoolingWeights(:))));
    end
    
    if (any(isnan(V1filterBank.sinPhasePoolingWeights)))
        fprintf('Sin pooling contains %d NANs\n', sum(isnan(V1filterBank.sinPhasePoolingWeights(:))));
    end
    
    V1filterBank.type = spatialPoolingKernelParams.type;
    V1filterBank.activationFunction = spatialPoolingKernelParams.activationFunction;
    V1filterBank.temporalPCAcoeffs = spatialPoolingKernelParams.temporalPCAcoeffs;
end

