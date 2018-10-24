function macularPigmentCalcs
    
    wavelength = 380:5:780;
    mDefault = Macular('wave', wavelength);

    densityAtEcc = Macular.eccDensity(0);
    mPeriphery = Macular('wave', wavelength, 'density', densityAtEcc);
    
    x = -2:0.1:2;
    xCenter = (numel(x)-1)/2
    y = x;
    [X,Y] = meshgrid(x,y);
    eccDegs = sqrt(X.^2+Y.^2);
    eccDegsVec = eccDegs(:);
    
    opticalImageRows = size(X,1);
    opticalImageCols = size(X,2);
    wavesNum = numel(mDefault.wave);
    opticalImageBoostFactor = zeros(opticalImageRows, opticalImageCols, wavesNum);

    for eccIndex = 1:numel(eccDegsVec)
        % Set the density
        mPeriphery.density = Macular.eccDensity(eccDegsVec(eccIndex));
        
        % Compute boost factor for optical image at this eccentricity
        [row,col] = ind2sub(size(X), eccIndex);
        opticalImageBoostFactor(row,col,:) = mPeriphery.transmittance./mDefault.transmittance;
    end
    
    for eccIndex = 1:size(X,1)
        
        ecc = x(eccIndex);
        % Set the density
        mPeriphery.density = Macular.eccDensity(ecc);
        
        figure(1);
        clf;
        subplot(2,2,1);
        plot(mDefault.wave,  mDefault.absorptance, 'k-'); hold on;
        plot(mPeriphery.wave,  mPeriphery.absorptance, 'rs-');
        plot(mPeriphery.wave, mPeriphery.transmittance, 'b-');
        ylabel(sprintf('absorptance\n(%% of quanta absorbed by MP))'));
        xlabel('wavelength (nm)');
    
        subplot(2,2,3);
        plot(mDefault.wave, mDefault.transmittance, 'ks'); hold on;
        plot(mDefault.wave, mPeriphery.transmittance, 'r'); hold on;
        legend({'transmitted photons (fovea)', sprintf('transmitted photons (%2.2f deg)', ecc)});
        ylabel('transmitted photons (%%)');
        xlabel('wavelength (nm)');
        
        subplot(2,2,4);
        plot(mDefault.wave, squeeze(opticalImageBoostFactor(xCenter,eccIndex,:)), 'ks');
        set(gca, 'YLim', [0.9 2.2]);
        ylabel('optical image boost factor');
        xlabel('wavelength (nm)');
        
        subplot(2,2,2);
        plot(mDefault.wave, mDefault.spectralDensity, 'k-');
        hold on;
        plot(mDefault.wave, mPeriphery.spectralDensity, 'rs-');
        title(sprintf('ecc = %f deg', ecc));
        ylabel('optical density');
        xlabel('wavelength (nm)');
        drawnow;
        pause()
    end
end

function Snodderly
    eccDegs = -5.0:0.01:5.0;
    [X,Y] = meshgrid(eccDegs, eccDegs);
    offset = 0.0;
    A = 0.64-offset;
    B = 0.42;
    eccDegsXY = sqrt(X.^2+Y.^2);
    mpDensitySnodderly = offset + A * 10.^(-B*eccDegsXY);
    
    mpDensityISETBio = 0.35 * 3.6028 ./ (eccDegsXY .^ 2 + 3.6028);
    
    mpDensity = mpDensityISETBio;
    
    mpDensityVec = mpDensity(:);

    
    fieldSizeDegsExamined = 1.0:0.1:20;
    for k = 1:numel(fieldSizeDegsExamined)
        fieldSizeDegs = fieldSizeDegsExamined(k);
%         idx = find(...
%             (abs(X(:)) < 0.5*fieldSizeDegs) & ...
%             (abs(Y(:)) < 0.5*fieldSizeDegs) );
        idx = find(eccDegsXY(:) <= 0.5*fieldSizeDegs);
        mpDensityField(k) = mean(mpDensityVec(idx));
    end
    
    figure(2); clf;
    subplot(1,2,1);
    plot(eccDegsXY(:), mpDensityVec, 'k.');
    set(gca, 'YLim', [0 0.7]);
    xlabel('eccentricity (degrees)');
    
    subplot(1,2,2);
    plot(fieldSizeDegsExamined, mpDensityField, 'k-');
    hold on;
    plot(fieldSizeDegsExamined, 0.485*exp(-fieldSizeDegsExamined/6), 'ro');
    xlabel('field size (degrees)');
end

