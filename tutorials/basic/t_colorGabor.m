function validationData = t_colorGabor(varargin)
% T_COLORGABOR  Calculate cone isomerizations for a static color Gabor modulation.
%     validationData = T_COLORGABOR(varargin)
% 
%     Illustrates the basic steps required to calculate cone isomerizations
%     for a static color Gabor modulation.
% 
%     Create a scene with a color gabor patch with color directions
%     specified as L, M, and S cone contrasts.  The scene will produce
%     a Gabor with these contrasts on a specified monitor.  Then passes the
%     scene through the optics and a cone mosaic and gets the isomerizations at
%     each cone.
% 
%     If parameters structure is not passed, the routine will use the defaults
%     provided by
%       responseParamsGenerate
%     That function and its subfunctions also documents what the relavant parameters are.
% 
%     The returned validation structure allows this routine to be called from a
%     validation script driven by the UnitTest toolbox.
% 
%     The tutorial produces output according to a scheme controlled by the
%     specified IBIOColorDetect rwObject.
% 
%     Optional key/value pairs
%      'rParams' - Value the is the rParams structure to use.  Default empty,
%         which then uses defaults produced by generation function.
%      'generatePlots' - true/false (default true).  Make plots?
%      'setRngSeed' - true/false (default true). When true, set the rng seed so noise is frozen.
%      'hexMosaic' - true/false (default false). Use a hexagonal mosaic, rather than a rectangular mosaic.
%      'wavelengths' - vector (default [400 10 700]). Start, delta, end wavelength sampling.
%        This parameter (including its default) overrides what is in
%        rParams, whether that is passed or obtained via
%        RESPONSEPARAMSGENERATE.
%
%     See also T_CONEISOMERIZATIONSMOVIE RESPONSEPARAMSGENERATE COLORSCENECREATE

%% Parse vargin for options passed here
p = inputParser;
p.addParameter('rParams',[],@isemptyorstruct);
p.addParameter('generatePlots',true,@islogical);
p.addParameter('setRngSeed',true,@islogical);
p.addParameter('hexMosaic',false,@islogical);
p.addParameter('wavelengths',[400 10 700],@isnumeric);
p.parse(varargin{:});
rParams = p.Results.rParams;

%% Clear
if (nargin == 0)
    ieInit; close all;
end

%% Fix random number generator so we can validate output exactly
if (p.Results.setRngSeed)
    rng(1);
end

%% Get the parameters we need
%
% Furnction responseParamsGenerate returns a hierarchical struct of
% parameters used by a number of tutorials and functions in this project.
% These provide reasonable defaults.
if (nargin < 1 | isempty(rParams))
    rParams = responseParamsGenerate;
    if (p.Results.hexMosaic)
        rParams.mosaicParams.conePacking = 'hex';
    end
end

%% Override wavelength parameters with specification
rParams.colorModulationParams.startWl = p.Results.wavelengths(1);
rParams.colorModulationParams.deltaWl = p.Results.wavelengths(2);
rParams.colorModulationParams.endWl = p.Results.wavelengths(3);

%% Set up the rw object for this program
rwObject = IBIOColorDetectReadWriteBasic;
theProgram = mfilename;
paramsList = {rParams.topLevelDirParams, rParams.mosaicParams, rParams.oiParams, rParams.spatialParams,  rParams.temporalParams,  rParams.backgroundParams, rParams.colorModulationParams};

%% Name the cone types for some printouts.
coneTypes = {'L' 'M' 'S'};
  
%% Create the isetbio scene for the gabor patch, using routine colorSceneCreate.
%
% This routine does a lot of colorimetry to produce a scene on a monitor
% that has the desired cone contrasts.
gaborScene = colorSceneCreate(rParams.spatialParams,rParams.backgroundParams,rParams.colorModulationParams,rParams.oiParams);

% Look at the scene image.  It is plausible for an L-M grating.  Remember that we
% are looking at the stimuli on a monitor different from the display file
% that we loaded, and thus the RGB values will not produce exactly the
% desired appearance.
if (p.Results.generatePlots)
    vcNewGraphWin; [~,h] = scenePlot(gaborScene,'radiance image no grid');
    rwObject.write('colorGaborScene',h,paramsList,theProgram,'Type','figure');
end

%% Create oi
gaborOI = oiCreate('wvf human');
gaborOI = oiSet(gaborOI,'h fov',rParams.spatialParams.fieldOfViewDegs);

% Set pupil diameter
focalLength = oiGet(gaborOI,'distance');
desiredFNumber = focalLength/(rParams.oiParams.pupilDiamMm/1000);
gaborOI  = oiSet(gaborOI ,'optics fnumber',desiredFNumber);
pupilDiamMmCheck = 1000*oiGet(gaborOI,'optics aperture diameter');
if (max(abs(pupilDiamMmCheck - rParams.oiParams.pupilDiamMm)) > 1e-8)
    error('Failed to set pupil diameter as expected');
end

%% Compute blurred optical image
%
% Turn of default off axis intensity falloff calculation first.
optics = oiGet(gaborOI,'optics');
optics = opticsSet(optics,'off axis method','skip');
gaborOI = oiSet(gaborOI,'optics',optics);
gaborOIBlur = oiCompute(gaborOI,gaborScene);

% Note how different the color appearance is than the scene.  This is
% because the OI incorprates the transmittance of the lens.  Down below we
% will turn that off as a check.
if (p.Results.generatePlots)
    vcNewGraphWin; [~,h] = oiPlot(gaborOIBlur,'irradiance image no grid');
    rwObject.write('colorGaborOpticalImageBlur',h,paramsList,theProgram,'Type','figure');
    clearvars('gaborOIBlur');
end

%% Turn off optics for current purpose of checking LMS contrast
% This involves replacing the OTF with a unity OTF, and recompute
optics = opticsSet(optics,'OTF',ones(size(opticsGet(optics,'OTF'))));
gaborOI = oiSet(gaborOI,'optics',optics);
gaborOI = oiCompute(gaborOI,gaborScene);

% Look at the OI
if (p.Results.generatePlots)
    vcNewGraphWin; [~,h] = oiPlot(gaborOI,'irradiance image no grid');
    rwObject.write('colorGaborOpticalImageNoBlur',h,paramsList,theProgram,'Type','figure');
end

% Just for fun, put OI into isetbio's interactive window
if (p.Results.generatePlots)
    vcAddAndSelectObject(gaborOI); oiWindow;
end

%% Verify that removing lens transmittance has expected effect
lens = oiGet(gaborOI,'lens');
lens.density = 0;
gaborOINoLens = oiSet(gaborOI,'lens',lens);
gaborOINoLens = oiCompute(gaborOINoLens,gaborScene);
if (p.Results.generatePlots)
    vcAddAndSelectObject(gaborOINoLens); oiWindow;
end

%% Create and get noise free sensor using coneMosaic obj
% Create a coneMosaic object here. When setting the fov, if only one value
% is specified, it will automatically make a square cone mosaic.
%
% You can generate either a hexagonal or a rectangular mosaic
if (strcmp(rParams.mosaicParams.conePacking, 'hex')) && ~any(isnan(rParams.mosaicParams.fieldOfViewDegs))
    % HEX mosaic
    resamplingFactor = 3;
    centerInMM = [0.0 0.0];                    % mosaic eccentricity
    spatiallyVaryingConeDensity = true;        % constant spatial density (at the mosaic's eccentricity)
    gaborConeMosaic = coneMosaicHex(resamplingFactor, spatiallyVaryingConeDensity, ...
        'center', centerInMM*1e-3, ...
        'spatialDensity', [0 rParams.mosaicParams.LMSRatio] ...
        );
    
    rParams.spatialParams.fieldOfViewDegs = 1;
    gaborConeMosaic.setSizeToFOVForHexMosaic(rParams.spatialParams.fieldOfViewDegs);
    gaborConeMosaic.visualizeGrid();
else
    % RECT mosaic
    gaborConeMosaic = coneMosaic;
    gaborConeMosaic.setSizeToFOV(rParams.spatialParams.fieldOfViewDegs);
end

% This used to be the default integration time, and we reset here to keep
% the vaidation happy.
gaborConeMosaic.integrationTime = 0.05;

% There is also an option of whether the cone current should be calculated
% in the compute function. If set to true, it uses an os object inside the
% coneMosaic object. The default is the linearOS.  Here we don't need that.
gaborConeMosaic.noiseFlag = 'none';
isomerizations = gaborConeMosaic.compute(gaborOI,'currentFlag',false);

%% Take a look at the mosaic responses in the window
if (p.Results.generatePlots)
    gaborConeMosaic.window;
end

% And must make a plot in a figure
if (p.Results.generatePlots)
    vcNewGraphWin; [~,h] = gaborConeMosaic.plot('cone mosaic');
    rwObject.write('colorGaborMosaic',h,paramsList,theProgram,'Type','figure');
    vcNewGraphWin; [~,h] = gaborConeMosaic.plot('mean absorptions');
    rwObject.write('colorGaborIsomerizations',h,paramsList,theProgram,'Type','figure');
end

%% Get min max for LMS cone isomerizations
% Extract the min and max absorptions in a loop. Since we are
% extracting only L, M, or S absorptions at each iteration, we get a vector
% so one call to max/min will suffice.
%
% These are close enough to the desired values that it seems OK, although
% why they aren't exactly the desired values is a little mysterious.  It is
% possible that the oi/coneMosaic object leads to slightly different cone
% fundamentals, or that the monitor quantization leads to the small
% deviatoins.  Someone energetic could track this down.
conePattern = gaborConeMosaic.pattern;
for ii = 2:4
    maxIsomerizations(ii) = max(isomerizations(conePattern==ii));
    minIsomerizations(ii) = min(isomerizations(conePattern==ii));
    contrasts(ii) = ...
        (maxIsomerizations(ii)-minIsomerizations(ii))/(maxIsomerizations(ii)+minIsomerizations(ii));
    
    fprintf('%s cone isomerzations\n\tMax: %d \n\tMin: %d\n',coneTypes{ii-1},maxIsomerizations(ii),minIsomerizations(ii));
    fprintf('\tAbsolute contrast: %04.3f\n',contrasts(ii));
end

%% Send back some validation data if requested
if (nargout > 0)
    validationData.maxIsomerizations = maxIsomerizations;
    validationData.minIsomerizations = minIsomerizations;
    validationData.contrasts = contrasts;
end

