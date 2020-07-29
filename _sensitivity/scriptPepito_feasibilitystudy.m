clear all;close all;
savePath = '/home/omartin/Projects/Turbulence Profiling/PEPITO/Simulations/';
%% Define the system
% Define an atmosphere
r0 = 0.1;
L0 = 25;
alt= 8e3;
w0 = 1;
wS = 15;
wD = pi/4;
atm = atmosphere(photometry.V0,r0,L0,'altitude',alt,'fractionnalR0',w0,...
    'windSpeed',wS,'windDirection',wD);

% Define a telescope
D   = 0.5;
fov = 2;
nRes= 100;
Fe  = 1000;
cobs= 0.39;
tel = telescope(D,'fieldOfViewInArcMin',fov,'resolution',nRes,...
    'samplingTime',1/Fe,'obstructionRatio',cobs);

% Define sources
n       = 9;
xx      = linspace(-fov*60/sqrt(2)/2,fov*60/sqrt(2)/2,n);               %
[xS,yS] = meshgrid(xx);
[zS,aS] = psfrTools.arcsec2polar(xS,yS);
ngs     = source('wavelength',photometry.V0,'zenith',zS,'azimuth',aS);
nSrc    = numel(ngs);
%% Define an imaging camera
expTime = 10e-3;
cam = imager('fieldStopSize',nRes/2,'nyquistSampling',1);
ngs = ngs.*tel*cam;
cam.referenceFrame = cam.frame;
cam.clockRate    = 1;
cam.exposureTime = expTime*Fe;
cam.startDelay   = 0;
flush(cam);
resCam = cam.resolution(1);
%% LOOP INSTANTIATION
tel = tel + atm;

nIter    = 10000;
nExp     = nIter/(Fe*expTime);
shortImg = zeros(resCam,resCam,nSrc,nExp);

%% RUN THE LOOP
reset(tel);
hwait = waitbar(0,'Loop is being closed...');
dt    = 0;
count = 0;
kCam  = 0;

for kIter=1:nIter
    tic;
    +tel;
    ngs = ngs.*tel*cam;
    
    if mod(kIter,Fe*expTime) == 0
        kCam = kCam + 1;
        shortImg(:,:,:,kCam) = reshape(cam.frame,resCam,resCam,nSrc);
    end
    
    % progress bar
    waitbar(kIter/nIter);    
    dt  = dt + toc();    
    if (kIter>1) && (kIter<nIter)
        fprintf(1, repmat('\b',1,count)); %delete line before
        count = fprintf('Remaining simulation time: %0.5g s',dt*(nIter-kIter)/kIter);
    elseif kIter==nIter
        fprintf(1, repmat('\b',1,count)); %delete line before
        count = fprintf('Remaining simulation time: %0.5g s\n',dt*(nIter-kIter)/kIter);
    end
end
close(hwait);

fitswrite(shortImg,[savePath,'simulatedShortExposureImages.fits']);

%% PROCESSING

% Load data
if ~exist('shortImg','var')
  shortImg = fitsread([savePath,'simulatedShortExposureImages.fits']);
end

resCam = size(shortImg,1);
nSrc   = size(shortImg,3);
nExp   = size(shortImg,4);

% Recenter the on-axis PSF
psfOn = squeeze(shortImg(:,:,floor(nSrc/2+1),:));
xOn   = zeros(1,nExp);
yOn   = zeros(1,nExp);

for k=1:nExp
    % Select the exposure
    tmp = squeeze(psfOn(:,:,k));
    % Get the barycenter
    [xOn(:,k),yOn(:,k)] = psfrTools.getBarycenter(tmp);
end

%% Apply the recentering to all PSF
shortImg_c = zeros(resCam,resCam,nSrc,nExp);
for k=1:nExp
    for j=1:nSrc
        % Select the exposure
        tmp = squeeze(shortImg(:,:,j,k));
        % Get the barycenter
        shortImg_c(:,:,j,k) = psfrTools.translateImage(tmp,[resCam/2-yOn(:,k),resCam/2-xOn(:,k)]);
    end
end

fitswrite(shortImg_c,[savePath,'simulatedShortExposureImages_recentered.fits']);

%% Measure the FWHM countour in long exposure mode
longImg = mean(shortImg_c,4);
fitswrite(longImg,[savePath,'simulatedLongExposureImages_recentered.fits']);

% PSF grid
psInMas = constants.radian2mas*ngs(1).wavelength/D/4;

% FWHM grid
for i=1:nSrc
    stat(i) = psfStats(longImg(:,:,i),tel.pupil,D,zS(i),aS(i),ngs(i).wavelength,psInMas,expTime);
end

% PSF grid
psfrTools.displayContours(stat,16);
