%% MANAGE workspaces
close all; % close all figures
clear all; % clear the workspace

path_oomao = '/home/omartin/Projects/SIMULATIONS/OOMAO/mfiles/'; % the oomao path
path_workspace = '/home/omartin/Projects/PROFILING/PEPITO/CODES/_processData/'; % the simulation folder path
path_pepito= '/home/omartin/Projects/PROFILING/PEPITO/pepito/'; % the simulation folder path
addpath(genpath(path_oomao),path_workspace,path_pepito);

%% GET DATA ID
path_data = '/run/media/omartin/HDD_OBM/PEPITO_DATA/'; % the simulation folder path
folder_data = dir(path_data);
dataID = {folder_data(contains({folder_data.name},'.fits')).name};
nObj = numel(dataID);
%% STARS positions

xStars = [122 158 449 524 54 449];
yStars = [703 708 464 325 273 46];
nStars = 4;%numel(xStars);
nCrop = 42;

%% HARDWARE CONFIGURATION
D = 0.5;
obs = 0.47;
wvl = 525e-9;
psInMas = 495.6;
Samp = constants.radian2mas*wvl/D/2/psInMas;
dk  = ceil(nCrop/Samp);
tel = telescope(D,'resolution',dk,'obstructionRatio',obs);


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
    %cube_LE(:,:,j) = tools.processImage(im_LE(idxj,idyj),0,1,0,Samp,...
                    %'fovInPixel',nCrop,'masking',false,'rebin',4,'tfccd',false,...
                    %'thresholding',-Inf);
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
        
