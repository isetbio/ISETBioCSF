function theOI = oiWithAOoptics(theOI)
    
    % Get the human focal length
    humanFocalLengthInMeters = oiGet(theOI,'focal plane distance');
    fNumber = oiGet(theOI, 'fNumber');
    wave = oiGet(theOI, 'optics wave');
    
    % Create diffraction limited oi
    theOI = oiCreate('diffraction limited');
    
    % Pull out the optics structure
    optics = oiGet(theOI,'optics');
    
    % Set the fLength and the fNumber to those of human
    optics = opticsSet(optics,'flength',humanFocalLengthInMeters); 
    optics = opticsSet(optics,'fnumber',fNumber);
    theOI = oiSet(theOI,'optics',optics);
 
    % Add human lens
    theOI = oiSet(theOI, 'lens', Lens('wave', wave));
end
