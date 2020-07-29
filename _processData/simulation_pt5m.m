clear all;close all;clc;

path_oomao = '/home/omartin/Projects/SIMULATIONS/OOMAO/mfiles/'; % the oomao path
path_workspace = '/home/omartin/Projects/PEPITO/'; % the simulation folder path
path_res= '/home/omartin/Projects/PEPITO/Results/'; % the simulation folder path

addpath(genpath(path_oomao),genpath(path_workspace));

%% Define the system
% Define the sources
psInMas = 495.6;
wvl     = 500e-9;

xStars_0 = ([122 158 449 524] - 608/2) * psInMas/(1e3 * 206264.8)/8; %in rad
yStars_0 = ([703 708 464 325] - 968/2) * psInMas/(1e3 * 206264.8)/8;
[aS,zS] = cart2pol(xStars_0,yStars_0);
ngs     = source('wavelength',photometry.V0,'zenith',zS,'azimuth',aS);
nSrc    = numel(ngs);
fov     = 2*max(hypot(xStars_0,yStars_0))*180*60/pi; %in arcmin


% Define the atmosphere
r0 = 0.1;
L0 = 16;
alt= [0,10e3];
w0 = [.7,.3];
wS = [5,20];
wD = [0,pi/4];

atm = atmosphere(photometry.V0,r0,L0,'altitude',alt,'fractionnalR0',w0,...
    'windSpeed',wS,'windDirection',wD);

% Define the telescope
D    = 0.5;
samp = 1e3*180*3600/pi*wvl/D/psInMas;
nPSF = 26;
k_   = fix(ceil(2.0/samp));
fovCam = nPSF * k_;
nRes = fovCam * 2;
Fe   = 1/0.04;
cobs = 0.47;
tel  = telescope(D,'fieldOfViewInArcMin',fov,'resolution',nRes,...
    'samplingTime',1/Fe,'obstructionRatio',cobs);


% Define the imaging camera
% fovCam is in wvl/D units:
% fovCam = nPSF * psInMas * (pi/(1e3 * 180 * 3600)) * (D/wvl)
% fovCam = nPSF * samp
% it's an integer -> fovCam = nPSF * k_/2 = nRes/2

cam = imager('nyquistSampling',1,'fieldStopSize',fovCam);



%% LOOP INSTANTIATION

% camera calibration
ngs = ngs.*tel*cam;
cam.referenceFrame = cam.frame;
cam.startDelay   = 0;
flush(cam);
nExp    = 50;
imSE = zeros(nPSF,nPSF,nSrc,nExp);

% phase propagation
tel = tel + atm;

%% RUN THE LOOP
reset(tel);

for kIter = 1 : nExp
    kIter
    +tel;
    ngs = ngs .* tel * cam;    
    imSE(:,:,:,kIter) = reshape(tools.binning(tools.crop(cam.frame,size(cam.frame)/2),[nPSF,nPSF*nSrc]),nPSF,nPSF,nSrc);        
end
psInMas = fovCam*wvl/D*constants.radian2mas/2/nPSF;

%% PEPITO
fitCn2 = false;
pep = pepito('parFile_pt5m','fitCn2',fitCn2);

% seeing measurements
psf_LE = squeeze(sum(imSE,4));
pep = pep.measureSeeing(psf_LE,'method','PSF','fitVib',false,'fitL0',false,'fitBg',false);
            
figure;
clf;
semilogy(pep.ydata);
hold on;
semilogy(pep.mod(pep.xfinal,pep.xdata));
