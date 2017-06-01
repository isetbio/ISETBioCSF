function theOI = updateOIWithCustomOptics(theOI, pupilDiameterMM, varargin)
% 
% Parse input
p = inputParser;
p.addRequired('theOI',@isstruct);
p.addRequired('pupilDiameterMM', @isnumeric);
p.addParameter('opticsModel', '', @ischar);
p.parse(theOI,pupilDiameterMM, varargin{:});

opticsModel = p.Results.opticsModel;
if (isempty(opticsModel))
    printf('Using default optics model\n');
    return;
elseif ismember(opticsModel, {'WvfHumanMeanOTF'})
else
    error('Unknown optics model: ''%s''.', opticsModel)
end

% Load the 4.5 mm pupil diameter Thibos data set
measPupilDiameterMM = 4.5;
[Zcoeffs_SampleMean, Zcoeffs_S] = wvfLoadThibosVirtualEyes(measPupilDiameterMM);
 

subjectsNum  = 100;
% Generate Z-coeffs for subjectsNum
Zcoeffs = ieMvnrnd(Zcoeffs_SampleMean, Zcoeffs_S, subjectsNum)'; 
    
phaseMethod = 'subject mean';

% Pull out the optics structure
optics = oiGet(theOI,'optics');

% Get the wavelengths
wavelengthsListToCompute = opticsGet(optics,'wave');

% Get the OTF sf support in c/deg
% xSfGridCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
% ySfGridCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
        
% sfValuesCyclesMm = opticsGet(optics,'otf support','mm');
% uMPerMm = 1000;
% uMPerDeg = 300;
% [xSfGridCyclesMm, ySfGridCyclesMm] = meshgrid(sfValuesCyclesMm{1},sfValuesCyclesMm{2});
% xSfGridCyclesDeg = uMPerDeg*xSfGridCyclesMm/uMPerMm;
% ySfGridCyclesDeg = uMPerDeg*ySfGridCyclesMm/uMPerMm;

% Get the psf spatial support from the spatial frequency support
% centerPosition = floor(length(xSfGridCyclesDeg)/2)+1;
% [xGridMinutes,yGridMinutes] = SfGridCyclesDegToPositionGridMinutes(xSfGridCyclesDeg,ySfGridCyclesDeg);
% position1DMinutes = xGridMinutes(centerPosition,:);

for waveIndex = 1:numel(wavelengthsListToCompute)
   wavelength = wavelengthsListToCompute(waveIndex);
    
   customMTFdata = computeCustomMTF(Zcoeffs, ...
       pupilDiameterMM, wavelength, subjectsNum, phaseMethod);
    
   if (waveIndex == 1)
        [X,Y] = meshgrid(customMTFdata.xSfGridCyclesDeg,customMTFdata.ySfGridCyclesDeg);
        [xGridMinutes,yGridMinutes] = SfGridCyclesDegToPositionGridMinutes(X,Y);
   end
   
   [a,b,thePsf] = OtfToPsf(X,Y,customMTFdata.otf);
   
   if (any(a ~= xGridMinutes | b ~= yGridMinutes))
      error('Internal coordinate system transformation error');
   end
        
   % Make sure psf has unit volume
   thePsf = thePsf/sum(thePsf(:));

   % Compute the OTF from the PSF
   [~,~,theOtfCentered] = PsfToOtf(xGridMinutes,yGridMinutes,thePsf);
   theOtf = ifftshift(theOtfCentered);
   if (waveIndex == 1)
       otfData = zeros(size(theOtf,1), size(theOtf,2), numel(wavelengthsListToCompute));
   end
   
   otfData(:,:,waveIndex) = theOtf;
end % waveIndex

% Update optics with new OTF data
optics = opticsSet(optics,'otf data',otfData);

% Stick optics into oi
theOI = oiSet(theOI,'optics',optics);

end

function customMTFdata = computeCustomMTF(Zcoeffs, pupilDiameterMM, wavelength, subjectsNum, phaseMethod)

    for subjectIndex = 1:subjectsNum
        % Make WVF for this subject
        wvfSubject = makeWVF(Zcoeffs(:,subjectIndex), wavelength, pupilDiameterMM, sprintf('subject-%d', subjectIndex));
        % Get PSF at selected wavelength
        psfSingleSubject = wvfGet(wvfSubject, 'psf', wavelength);
        optics = oiGet(wvf2oi(wvfSubject),'optics'); 
        % Get OTF at selected wavelength
        otf = opticsGet(optics, 'otf');
        sfx = opticsGet(optics, 'otf fx', 'cyclesperdeg');
        sfy = opticsGet(optics, 'otf fy', 'cyclesperdeg');
        otf = otfWithZeroCenteredPSF(otf, psfSingleSubject);
       
        if (subjectIndex == 1)
            otfSubjects = zeros(subjectsNum, size(otf,1), size(otf,2));
            customMTFdata.xSfGridCyclesDeg = sfx;
            customMTFdata.ySfGridCyclesDeg = sfy;
        end
        otfSubjects(subjectIndex,:,:) = otf;
    end % subjectIndex
    
    % Compute average MTF (mag and phase) across subjects
    otfMean = squeeze(mean(otfSubjects,1));
    otfMeanMag = fftshift(abs(otfMean));
    otfMeanPhase = fftshift(angle(otfMean));
    
    if strcmp(phaseMethod, 'subject mean')
       phaseToUse = otfMeanPhase;
    elseif strcmp(phaseMethod, 'zero')
       phaseToUse = 0*otfMeanPhase;
    end
        
    % Make radially symmetric magnitude
    centerPosition = floor(size(otfMeanMag,1)/2) + 1;
    otfSlice = squeeze(otfMeanMag(centerPosition, :));
    radiusMatrix = MakeRadiusMat(length(otfSlice),length(otfSlice),centerPosition);
    otfMeanMag = interp1(0:floor(length(otfSlice)/2),otfSlice(centerPosition:end),radiusMatrix,'linear',0);

    customMTFdata.otf = otfMeanMag .* exp(1i*phaseToUse);
end

function otf = otfWithZeroCenteredPSF(otf, psf)
    % Compute center of mass
    centerOfMass = computeCenterOfMass(psf);
   
    % Compute translation vector to bring center of mass at 0,0
    centerPosition = floor(size(psf,1)/2) + 1 * [1 1];
    translationVector = (centerPosition-centerOfMass);
    otf = shiftInFTplane(otf, translationVector);
end

function centerOfMass = computeCenterOfMass(psf0)
% Compute center of mass
    [rc,cc] = ndgrid(1:size(psf0,1),1:size(psf0,2));
    Mt = sum(psf0(:));
    centerOfMassY = sum(psf0(:) .* rc(:)) / Mt;
    centerOfMassX = sum(psf0(:) .* cc(:)) / Mt;
    centerOfMass = [centerOfMassX centerOfMassY];
end

function theWVF = makeWVF(zcoeffs, wavelengthsToCompute, pupilDiameterMM, name)
    theWVF = wvfCreate(...
                'wave',wavelengthsToCompute,...
                'zcoeffs',zcoeffs,...
                'name', name);
    theWVF = wvfSet(theWVF,'calc pupil size',pupilDiameterMM);
    theWVF = wvfComputePSF(theWVF);
end

function otf = shiftInFTplane(otf, translationVector)
    
    % The size of the matrix.
    [N, M] = size(otf);

    % The mathsy bit. The floors take care of odd-length signals.
    x_shift = exp(-1i * 2 * pi * translationVector(2) * [0:floor(N/2)-1 floor(-N/2):-1]' / N);
    y_shift = exp(-1i * 2 * pi * translationVector(1) * [0:floor(M/2)-1 floor(-M/2):-1] / M);

    % Force conjugate symmetry. Otherwise this frequency component has no
    % corresponding negative frequency to cancel out its imaginary part.
    if mod(N, 2) == 0
        x_shift(N/2+1) = real(x_shift(N/2+1));
    end 
    if mod(M, 2) == 0
        y_shift(M/2+1) = real(y_shift(M/2+1));
    end

    otf = otf .* (x_shift * y_shift);
end

