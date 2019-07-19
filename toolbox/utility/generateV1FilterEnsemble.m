function [V1filterEnsemble, hFig] = generateV1FilterEnsemble(spatialParams, mosaicParams, topLevelDirParams, visualizeSpatialScheme, thresholdParams, paramsList)

    originalStimulusFOV = spatialParams.fieldOfViewDegs;
    spatialParams.fieldOfViewDegs = max([spatialParams.fieldOfViewDegs mosaicParams.fieldOfViewDegs]);
    
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
        coneLocsInMeters = theMosaic.coneLocsHexGridAlignedWithSerializedConeMosaicResponse();  % OLD, incorrect way: coneLocsHexGrid; 
    end
    coneLocsDegs(:,1) = coneLocsInMeters(:,1) / theMosaic.width * theMosaic.fov(1);
    coneLocsDegs(:,2) = coneLocsInMeters(:,2) / theMosaic.height * theMosaic.fov(2);
       
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
    
    % modify some params for small V1 RFs
    spatialParams.windowType = 'Gaussian';
    
    cyclesPerRFs = thresholdParams.spatialPoolingKernelParams.cyclesPerRFs;
    orientationRFs = thresholdParams.spatialPoolingKernelParams.orientations;
    ensemblePositions = thresholdParams.spatialPoolingKernelParams.spatialPositionsNum;
    pixelsPerDeg = spatialParams.row/spatialParams.fieldOfViewDegs;
    ensembleSampleSpacingPixels = thresholdParams.spatialPoolingKernelParams.spatialPositionOffsetDegs * pixelsPerDeg;
    
    unitIndex = 0;
    ft2Dindex = 0;
    
    for bandwidthIndex = 1:numel(cyclesPerRFs)
        cyclesPerRF = cyclesPerRFs(bandwidthIndex);
        spatialParams.gaussianFWHMDegs = cyclesPerRF/spatialParams.cyclesPerDegree;
        
        v1Unit.cyclesPerRF = cyclesPerRF;
        v1Unit.bandwidthIndex = bandwidthIndex;
         
        v1Unit.rowsNum = 2*ensemblePositions+1;
        v1Unit.colsNum = 2*ensemblePositions+1;
        
        for orientationIndex = 1:numel(orientationRFs)
            orientationRF = orientationRFs(orientationIndex);
            v1Unit.orientationIndex = orientationIndex;
            v1Unit.orientationRF = orientationRF;
            spatialParams.ang = orientationRF/180*pi;
            
            ft2Dindex = ft2Dindex + 1;
            v1Unit.ft2Dindex = ft2Dindex;
            
            for rowOffset = -ensemblePositions:ensemblePositions
            for colOffset = -ensemblePositions:ensemblePositions
            
                spatialParams.center = [colOffset rowOffset]*ensembleSampleSpacingPixels;
                v1Unit.rowColPosition = [rowOffset colOffset];
                v1Unit.spatialPosition = spatialParams.center/(spatialParams.row/2) * originalStimulusFOV/2;
            
                % re-center RF to nearest cone
                %[~, nearestConeIndex] = min(sqrt(sum((bsxfun(@minus, coneLocsDegs, v1Unit.spatialPosition)).^2,2)));
                %v1Unit.spatialPosition = coneLocsDegs(nearestConeIndex,:);
                %spatialParams.center = v1Unit.spatialPosition * (spatialParams.row/2)/(spatialParams.fieldOfViewDegs/2);
            
                spatialParams.ph = 0;
                cosPhaseFilter = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, 1.0));
                v1Unit.cosPhasePoolingProfile = cosPhaseFilter-1;
                spatialParams.ph = pi/2;
                sinPhaseFilter = imageHarmonic(imageHarmonicParamsFromGaborParams(spatialParams, 1.0));
                v1Unit.sinPhasePoolingProfile = sinPhaseFilter-1;
              
                % Compute energy envelope
                RFprofile = (v1Unit.cosPhasePoolingProfile).^2 + (v1Unit.sinPhasePoolingProfile).^2;
                v1Unit.RFprofile = RFprofile / max(RFprofile(:));
            
                %if (rowOffset == -ensemblePositions) && (colOffset == -ensemblePositions)
                    xaxisDegs = (0:(size(RFprofile,2)-1))/size(RFprofile,2) * spatialParams.fieldOfViewDegs;
                    xaxisDegs = xaxisDegs - mean(xaxisDegs);
                    yaxisDegs = (0:(size(RFprofile,1)-1))/size(RFprofile,1) * spatialParams.fieldOfViewDegs;
                    yaxisDegs = yaxisDegs - mean(yaxisDegs);
                    [X,Y] = meshgrid(xaxisDegs+v1Unit.spatialPosition(1), yaxisDegs+v1Unit.spatialPosition(2));
                    rfCoordsDegs = [X(:) Y(:)];
                    % Find nearest cone location
                    [~, idx] = pdist2(rfCoordsDegs, coneLocsDegs, 'euclidean', 'Smallest', 1);
                %end
            
                % Sample according to cone locations
                v1Unit.cosPhasePoolingWeights = v1Unit.cosPhasePoolingProfile(idx);
                v1Unit.sinPhasePoolingWeights = v1Unit.sinPhasePoolingProfile(idx);
                
                % Adjust weights by the inverse of the coneDensity
                if (thresholdParams.spatialPoolingKernelParams.adjustForConeDensity)
                    v1Unit.cosPhasePoolingWeights = v1Unit.cosPhasePoolingWeights ./ coneDensity;
                    v1Unit.sinPhasePoolingWeights = v1Unit.sinPhasePoolingWeights ./ coneDensity;
                end
    
                v1Unit.envelopePoolingWeights = (v1Unit.cosPhasePoolingWeights).^2 + (v1Unit.sinPhasePoolingWeights).^2;

                v1Unit.type = thresholdParams.spatialPoolingKernelParams.type;
                v1Unit.activationFunction = thresholdParams.spatialPoolingKernelParams.activationFunction;
                v1Unit.temporalPCAcoeffs = thresholdParams.spatialPoolingKernelParams.temporalPCAcoeffs;

                unitIndex = unitIndex + 1;
                V1filterEnsemble{unitIndex} = v1Unit;
            end
            end
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
        paramsList{numel(paramsList)+1} = thresholdParams;
        if (numel(hFig) == 1)
            fileName = sprintf('V1poolingKernel');
            rwObject.write(fileName, data, paramsList, theProgram, ...
            'type', 'NicePlotExportPDF', 'FigureHandle', hFig, 'FigureType', 'pdf');
        else
            for figIndex = 1:numel(hFig)
                theFigHandle = hFig(figIndex);
                fileName = sprintf('V1poolingEnsemble_%d', figIndex);
                rwObject.write(fileName, data, paramsList, theProgram, ...
                    'type', 'NicePlotExportPNG', 'FigureHandle', theFigHandle , 'FigureType', 'png');
            end
        end
    end % visualizeSpatialScheme
end
