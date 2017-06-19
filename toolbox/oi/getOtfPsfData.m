function [otf, otf_fx, otf_fy, psf, psf_x, psf_y] = getOtfPsfData(theOI, wavelengthsToCompute)
    theOptics = oiGet(theOI, 'optics');
    otf = abs(fftshift(opticsGet(theOptics,'otf data', wavelengthsToCompute)));
    otf_fx = opticsGet(theOptics, 'otf fx', 'cyclesperdeg');
    otf_fy = opticsGet(theOptics, 'otf fx', 'cyclesperdeg');
    
    psf = opticsGet(theOptics,'psf data', wavelengthsToCompute);
    %psfSupport = opticsGet(theOptics,'psf support', 'deg');
    [X,Y] = meshgrid(otf_fx,otf_fy);
    [xGridMinutes,yGridMinutes] = SfGridCyclesDegToPositionGridMinutes(X,Y);
    psf_x = xGridMinutes(1,:);
    psf_y = psf_x;
end
