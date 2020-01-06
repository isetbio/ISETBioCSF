function makeZones()
    eccentricitiesInDegs = 0:0.05:20;
    eccentricitiesInMicrons = eccentricitiesInDegs * 300;
    eccentricitiesInMeters = eccentricitiesInMicrons * 1e-6;
    angles = 0*eccentricitiesInMeters;
    
    [coneSpacingInMeters, aperture, density] = ...
        coneSizeReadData('eccentricity', eccentricitiesInMeters, ...
        'angle', angles);
    coneSpacingInMicrons1 = coneSpacingInMeters' * 1e6;
    
    [coneSpacingInMeters, aperture, density] = ...
        coneSizeReadData('eccentricity', eccentricitiesInMeters, ...
        'angle', angles+90);
    coneSpacingInMicrons2 = coneSpacingInMeters' * 1e6;
    
    [coneSpacingInMeters, aperture, density] = ...
        coneSizeReadData('eccentricity', eccentricitiesInMeters, ...
        'angle', angles+180);
    coneSpacingInMicrons3 = coneSpacingInMeters' * 1e6;
    
    [coneSpacingInMeters, aperture, density] = ...
        coneSizeReadData('eccentricity', eccentricitiesInMeters, ...
        'angle', angles+270);
    coneSpacingInMicrons4 = coneSpacingInMeters' * 1e6;
    
    averageSpacing = (coneSpacingInMicrons1+coneSpacingInMicrons2+coneSpacingInMicrons3+coneSpacingInMicrons4)/4;
    minSpacing = min(averageSpacing);
    
    for spacingStep = 1:15
        idx = find(averageSpacing>=minSpacing & averageSpacing < minSpacing+1);
        ecc = eccentricitiesInMicrons(idx);
        eccRange{spacingStep} = [min(ecc)  max(ecc)];
        fprintf('Spacing step%d, ecc:[%2.3f - %2.3f]\n', spacingStep, min(ecc),  max(ecc));
        minSpacing = minSpacing + 1;
    end
    
    figure();
    plot(eccentricitiesInMicrons, coneSpacingInMicrons1, 'm.'); hold on
    plot(eccentricitiesInMicrons, coneSpacingInMicrons2, 'r.');
    plot(eccentricitiesInMicrons, coneSpacingInMicrons3, 'b.');
    plot(eccentricitiesInMicrons, coneSpacingInMicrons4, 'g.');
    plot(eccentricitiesInMicrons, averageSpacing, 'k-', 'LineWidth', 2);
    for spacingStep = 1:15
        range = eccRange{spacingStep}
        if (~isempty(range))
            plot(range(1)*[1 1], [0 20], 'k-');
        end
    end
    grid on;
    set(gca, 'YTick', 0:1:20)
    xlabel('ecc (microns)');
    ylabel('spacing (microns)');
end
