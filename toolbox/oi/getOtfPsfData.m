% Method to get the OTF and the PSF for a particular wavelength (as well as the microns/degree factor) from an oi
function [otf, otf_fxCyclesPerDeg, otf_fyCyclesPerDeg, psf, psf_xMinutes, psf_yMinutes, micronsPerDegree] = getOtfPsfData(theOI, selectedWavelength)
    % Get OTF
    theOptics = oiGet(theOI, 'optics');
    otf = abs(fftshift(opticsGet(theOptics,'otf data', selectedWavelength)));
    otf_fxCyclesPerDeg = opticsGet(theOptics, 'otf fx', 'cyclesperdeg');
    otf_fyCyclesPerDeg = opticsGet(theOptics, 'otf fx', 'cyclesperdeg');
    [xSfGridCyclesDegGrid,ySfGridCyclesDegGrid] = meshgrid(otf_fxCyclesPerDeg, otf_fyCyclesPerDeg);
    
    % Get PSF and its support (in minutes) via OtfToPsf
    [xGridMinutes,yGridMinutes,psf] = OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,otf);
    psf_xMinutes = squeeze(xGridMinutes(1,:));
    psf_yMinutes = squeeze(yGridMinutes(:,1));
            
    % Get psf support in microns
    psfMicrons = opticsGet(theOptics,'psf support','um');
    psfXsupportMicrons = psfMicrons{1}(1,:);

    % Compute microns-per-degree for this OI
    psfXsupportDegrees = psf_xMinutes*60;
    micronsPerDegree = psfXsupportMicrons(1)/psfXsupportDegrees(1);
end
