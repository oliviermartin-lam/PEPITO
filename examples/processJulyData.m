%% MANAGE workspaces
close all; % close all figures
clear all; % clear the workspace
path_pepito = '/home/omartin/Projects/PEPITO/PEPITO/';
addpath(genpath(path_pepito));

%% GET DATA ID
path_data = '/run/media/omartin/HDD_OBM/PEPITO_DATA/'; % the data folder path
folder_data = dir(path_data);
dataID = {folder_data(contains({folder_data.name},'.fits')).name};
nObj = numel(dataID);

%% STARS positions
xStars = [122 158 449 524 54 449]; % Stars x position in pixels
yStars = [703 708 464 325 273 46]; % Stars y positions in pixels
nStars = 4;%numel(xStars);
nCrop = 43; 
% Note that there is no automatic detection yet, which is not
% pipeline-compliant

%% HARDWARE CONFIGURATION
D = 0.5; % telescope diameter in m
obs = 0.47; % telescope central obscuration
wvl = 525e-9; % filter central wavelength (to be verified)
psInMas = 495.6; % detector pixel scale (to be verified)
Samp = constants.radian2mas*wvl/D/2/psInMas; %PSF sampling in lambda/(2D) units
dk  = 2*round(nCrop/Samp/2); % telescope resolution:
%Note that for a PSF nyquist-sampled PSF (Samp=1), the pupil resolution must be twice lower than the PSF/OTF size as the OTF is the pupil auto-correlation

% Define Pupil
x = linspace(-D/2,D/2,dk);
[X,Y] = meshgrid(x);
P = hypot(X,Y) < D/2 & hypot(X,Y)>=obs*D/2;
tel = struct('D',D,'resolution',dk,'obs',obs,'pupil',P); % in the future, we may think about having the telescope structure defintion within pepito directly


%% GRAB DATA
% short exposure images
k = 7;
im_SE = fitsread([path_data,cell2mat(dataID(k))]);
% long exposure images
im_LE = sum(im_SE,3);
cube_LE = zeros(nCrop,nCrop,nStars);

for j=1:nStars
    idxj = xStars(j) - nCrop/2+1:xStars(j)+nCrop/2;
    idyj = yStars(j) - nCrop/2+1:yStars(j)+nCrop/2;
    cube_LE(:,:,j) = im_LE(idxj,idyj);    
end

[x0,y0,nSE] = size(im_SE);

fInit =zeros(1,nStars);
for k=1:nStars
    fInit(k) = sum(sum(cube_LE(:,:,k)));
end
fInit = fInit/sum(fInit(:));
    
%% RUN PEPITO

pep = pepito(im_SE,tel,psInMas,xStars(1:nStars),yStars(1:nStars),nCrop,3,'wvl',wvl);

%% 1. Seeing estimation

% Get long exposure images
pep.LEframes = pep.stackingFrames(pep.SEframes,3);
% Croping
pep.SEpsf = pepito.extractPSF(pep.SEframes,pep.xStars,pep.yStars,pep.nBox);
% Stacking
pep.LEpsf = pepito.stackingFrames(pep.SEpsf,3);
% Estimating the seeing
% pep = pep.measureSeeing(pep.LEpsf,'weighting',false);%,'vibInit',[]);
% 
% % Display
% ims = reshape(pep.im_sky,nCrop,nCrop*nStars);
% imf = reshape(pep.im_fit,nCrop,nCrop*nStars);
% imr = reshape(pep.im_fit   - pep.im_sky,nCrop,nCrop*nStars);
% h = figure;
% imagesc(asinh([ims;imf;imr]));

%2.1539    2.8800    6.5673    2.6006   21.2941
%% 2. Cn2h estimation

%1\ Get anisoplanatic PSF
% Measure tip-tilt
[xOn,yOn] = pepito.getBarycenter(pep.SEpsf);          
% Compensate tip-tilt
SEframes_c = pepito.recenterFrames(pep.SEpsf,yOn,xOn);
% Long-exposure
out        = pepito.stackingFrames(SEframes_c,3);

%%
kim = 1;kref=4;nB = 16;
a = pep.LEpsf(:,:,kim)/sum(sum(pep.LEpsf(:,:,kim)));
a = a(19-nB/2+1:19+nB/2,20-nB/2+1:20+nB/2);
b = squeeze(out(:,:,kim,kref));
b = b(19-nB/2+1:19+nB/2,20-nB/2+1:20+nB/2);
close all;figure,imagesc([a,b]);
set(gca,'XTick',[],'YTick',[]);
pbaspect([2,1,1]);
        
