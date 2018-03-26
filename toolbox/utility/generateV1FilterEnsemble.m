function [V1filterEnsemble, hFig] = generateV1FilterEnsemble(spatialParams, mosaicParams, topLevelDirParams, visualizeSpatialScheme, thresholdParams, paramsList)

    % Make the stimulus spatial modulation
    spatialPattern = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, 1.0));
    spatialModulation = spatialPattern-1;
    spatialModulation = spatialModulation/max(abs(spatialModulation(:)));
            
    % Load the mosaic
    coneParamsList = {topLevelDirParams, mosaicParams};
    theProgram = 't_coneCurrentEyeMovementsResponseInstances';
    rwObject = IBIOColorDetectReadWriteBasic;
    theMosaic = rwObject.read('coneMosaic', coneParamsList, theProgram, 'type', 'mat');
    
    % Get the cone locs in degrees
    if (strcmp(mosaicParams.conePacking,'rect'))
        coneLocsInMeters = theMosaic.coneLocs;
    else
        coneLocsInMeters = theMosaic.coneLocsHexGrid;
    end
    coneLocsDegs(:,1) = coneLocsInMeters(:,1) / theMosaic.width * theMosaic.fov(1);
    coneLocsDegs(:,2) = coneLocsInMeters(:,2) / theMosaic.height * theMosaic.fov(2);
       
    if (strcmp(mosaicParams.conePacking, 'hex'))
        % Find the density map around each cone
        eccInMeters = sqrt(sum(coneLocsInMeters.^2, 2));
        ang = atan2(squeeze(coneLocsInMeters(:,2)), squeeze(coneLocsInMeters(:,1)))/pi*180;
        [~, ~, coneDensity] = coneSizeReadData('eccentricity',eccInMeters(:),'angle',ang(:));
    else
        coneDensity = 1;
    end
    
    % modify some params for small V1 RFs
    spatialParams.windowType = 'Gaussian';
    cyclesPerRF = 1.5;
    spatialParams.gaussianFWHMDegs = cyclesPerRF/spatialParams.cyclesPerDegree;
    
    ensemblePositions = 4;
    ensembleSampleSpacing = round((spatialParams.row/2)/ensemblePositions);
    
    fprintf('Generating %d v1 pooling units\n', (2*ensemblePositions+1)^2);
    unitIndex = 0;
    for rowOffset = -ensemblePositions:ensemblePositions
        for colOffset = -ensemblePositions:ensemblePositions
            
            spatialParams.center = [colOffset rowOffset]*ensembleSampleSpacing;
            v1Unit.spatialPosition = spatialParams.center * spatialParams.fieldOfViewDegs/2;
            v1Unit.rowColPosition = [colOffset rowOffset];
            
            spatialParams.ph = 0;
            cosPhaseFilter = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, 1.0));
            v1Unit.cosPhasePoolingProfile = cosPhaseFilter-1;
            spatialParams.ph = pi/2;
            sinPhaseFilter = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, 1.0));
            v1Unit.sinPhasePoolingProfile = sinPhaseFilter-1;
              
            % Compute energy envelope
            RFprofile = v1Unit.cosPhasePoolingProfile.^2 + v1Unit.sinPhasePoolingProfile.^2;
            v1Unit.RFprofile = RFprofile / max(abs(RFprofile(:)));
            
            if (rowOffset == -ensemblePositions) && (colOffset == -ensemblePositions)
                xaxisDegs = (0:(size(RFprofile,2)-1))/size(RFprofile,2) * spatialParams.fieldOfViewDegs;
                xaxisDegs = xaxisDegs - mean(xaxisDegs);
                yaxisDegs = (0:(size(RFprofile,1)-1))/size(RFprofile,1) * spatialParams.fieldOfViewDegs;
                yaxisDegs = yaxisDegs - mean(yaxisDegs);
                [X,Y] = meshgrid(xaxisDegs, yaxisDegs);
                rfCoordsDegs = [X(:) Y(:)];
                [~, idx] = pdist2(rfCoordsDegs, coneLocsDegs, 'euclidean', 'Smallest', 1);
            end
            
            v1Unit.outlineX = spatialParams.gaussianFWHMDegs * cosd(0:5:360);
            v1Unit.outlineX = spatialParams.gaussianFWHMDegs * sind(0:5:360);
            
            if (1==2)
                figure(1);
                subplot(1,3,1)
                imagesc(xaxisDegs, yaxisDegs, v1Unit.sinPhasePoolingProfile);
                set(gca, 'CLim', [-1 1]);
                axis 'square'
                subplot(1,3,2)
                imagesc(xaxisDegs, yaxisDegs, v1Unit.cosPhasePoolingProfile);
                set(gca, 'CLim', [-1 1]);
                axis 'square'
                subplot(1,3,3)
                imagesc(xaxisDegs, yaxisDegs, v1Unit.RFprofile);
                set(gca, 'CLim', [-1 1]);
                axis 'square'
                colormap(gray);
                drawnow;
            end
            
            v1Unit.cosPhasePoolingWeights = v1Unit.cosPhasePoolingProfile(idx);
            v1Unit.sinPhasePoolingWeights = v1Unit.sinPhasePoolingProfile(idx);
            v1Unit.envelopePoolingWeights = v1Unit.RFprofile(idx);
    
             
            % Adjust weights by the inverse of the coneDensity
            if (thresholdParams.spatialPoolingKernelParams.adjustForConeDensity)
                v1Unit.cosPhasePoolingWeights = v1Unit.cosPhasePoolingWeights ./ coneDensity;
                v1Unit.sinPhasePoolingWeights = v1Unit.sinPhasePoolingWeights ./ coneDensity;
                v1Unit.envelopePoolingWeights = v1Unit.envelopePoolingWeights ./ coneDensity;
            end
    
            v1Unit.type = thresholdParams.spatialPoolingKernelParams.type;
            v1Unit.activationFunction = thresholdParams.spatialPoolingKernelParams.activationFunction;
            v1Unit.temporalPCAcoeffs = thresholdParams.spatialPoolingKernelParams.temporalPCAcoeffs;
    
            unitIndex = unitIndex + 1;
            V1filterEnsemble{unitIndex} = v1Unit;
        end
    end
    

    hFig = [];
    if (visualizeSpatialScheme)
        coneRadiusMicrons = mosaicParams.innerSegmentSizeMicrons/2;
        coneRadiusMicrons = 0.5*diameterForCircularApertureFromWidthForSquareAperture(theMosaic.pigment.width)*1e6;
        hFig = visualizeSpatialPoolingScheme(xaxisDegs, yaxisDegs, spatialModulation, ...
            thresholdParams.spatialPoolingKernelParams, V1filterEnsemble, coneLocsDegs, ...
            mosaicParams.fieldOfViewDegs, spatialParams.fieldOfViewDegs, coneRadiusMicrons);
    
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
