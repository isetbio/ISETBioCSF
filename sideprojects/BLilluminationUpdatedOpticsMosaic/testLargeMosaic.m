function testLargeMosaic

    ibioDataDir = '/Volumes/DropBoxDisk/Dropbox/Dropbox (Aguirre-Brainard Lab)/IBIO_data';
    %ibioDataDir = '/media/dropbox_disk/dropboxlab/IBIO_data';
    %ibioDataDir = '/Volumes/SamsungT3/Dropbox/AguirreBrainardLabsDropbox/IBIO_data';
    %ibioDataDir = '/Volumes/IthakasPassport/DropboxLab/AquirreBrainardLab/IBIO_data'

    sceneRootDir = fullfile(ibioDataDir,'BLIlluminationDiscrimination/SceneData');
    oiRootDir = fullfile(ibioDataDir, 'BLIlluminationDiscrimination/OpticalImageDataNoBlur');
    
    %opticsModel = 'ThibosDefaultSubject3MMPupil';
    opticsModel = 'DeltaFunction';
    
    sceneNames = {...
        'Constant_CorrectSize'...
      %  'NM1_CorrectSize' ...		
      %  'NM2_CorrectSize' ...		
      %  'Neutral_CorrectSize' ...
        };
    illuminationNames = { ...
        'BlueIllumination' ...
        'GreenIllumination' ...	
        'RedIllumination' ...		
        'Standard' ...		
        'YellowIllumination' ...
        };
    %illuminationNames = {illuminationNames{4}};
    
    sceneNamesNum = 1; % numel(sceneNames)
    illuminationsNum = 1; % numel(illuminationNames)
    for sceneIndex = 1:sceneNamesNum
        sceneName = sceneNames{sceneIndex};
        for illumIndex = 1:illuminationsNum
        
            illuminationName = illuminationNames{illumIndex};
        
            sceneDir = fullfile(sceneRootDir, sceneName, illuminationName);
            oiDir = fullfile(oiRootDir, sceneName, illuminationName);
        
            listings = dir(sprintf('%s/*.mat', sceneDir));
        
            readOneOI = false;
            
            for k = 1:numel(listings)
                if (contains(listings(k).name, '.mat'))  && (~readOneOI)
                    % Source filename
                    
                    oiFileName = fullfile(oiDir, listings(k).name);
                    load(oiFileName, 'oi');
                    fprintf('oisize: %f x %f degs\n', oiGet(oi,'hfov'), oiGet(oi, 'vfov'));
                    readOneOI = true;
                end
            end % k
        end % illumIndex
    end % sceneIndex
    
    resamplingFactor = 4;
    integrationTime = 50/1000;
    
    if (1==2)
        cEccBased = coneMosaicHex(resamplingFactor, ...
            'fovDegs', [oiGet(oi,'hfov') oiGet(oi, 'vfov')], ...
            'integrationTime', integrationTime, ...
            'eccBasedConeDensity', true, ...
            'eccBasedConeQuantalEfficiency', true, ...
            'maxGridAdjustmentIterations', 40);
        fprintf('Cones num: %d\n', numel(find(cEccBased.pattern > 1)));
        
        save('largeMosaic.mat', 'cEccBased');
    else
        load('largeMosaic.mat', 'cEccBased');
        cEccBased.eccBasedConeQuantalEfficiency = true;
        %save('largeMosaic.mat', 'cEccBased');
    end
    
    
    
    cEccBased.displayInfo
    cEccBased.integrationTime
    figure(111);
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);
    cEccBased.visualizeGrid('axesHandle', ax);
    set(ax, 'YLim', [-1 1]*cEccBased.micronsPerDegree*1e-6, 'XLim', [-12 0]*cEccBased.micronsPerDegree*1e-6);
    drawnow;
    
    pause
    
    tic
    
    figure(111); clf;
    ax = subplot('Position', [0.05, 0.05, 0.94, 0.94]);
    isomerizations = cEccBased.compute(oi);
    toc
    
    size(isomerizations);
    signalRange = [min(isomerizations(:)) max(isomerizations(:))]
    showColorBar = true;
    labelColorBarTicks = true;
    backgroundColor = [0 0 0 ];
    xRange = [];
    yRange = [];
    
    cEccBased.renderActivationMap(ax, isomerizations, ...
             'signalRange', signalRange, ...
             'visualizedConeAperture', 'geometricArea', ...
             'mapType', 'modulated disks', ...
             'showColorBar', showColorBar, ...
             'labelColorBarTicks', labelColorBarTicks, ...
             'xRange', xRange, ...
             'yRange', yRange, ...
             'showXLabel', false, ...
             'showYLabel', false, ...
             'colorMap', gray(1024), ...
             'backgroundColor', backgroundColor);
        ylabel(ax, '');
        axis(ax, 'ij');
        xlim(ax, xRange*1e-6);
        ylim(ax, yRange*1e-6);
    cEccBased.displayInfo();
end
