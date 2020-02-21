close all; % close all figures
clear all; % clear the workspace
path_pepito = '/home/omartin/Projects/PEPITO/PEPITO/';
addpath(genpath(path_pepito));

%% GET DATA ID
path_data = '/run/media/omartin/8C8C-B02C/PEPITO_DATA/'; % the data folder path
folder_data = dir(path_data);
dataID = {folder_data(contains({folder_data.name},'.fits')).name};
nObj = numel(dataID);
%% STARS positions
xStars = [122 158 449 524 54 449]; % Stars x position in pixels
yStars = [703 708 464 325 273 46]; % Stars y positions in pixels
nStars = 4;%numel(xStars);
nCrop = 43; 
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
tel = struct('D',D,'resolution',dk,'obs',obs,'pupil',P); 

%% GRAB DATA
% short exposure images
k = 7;
im_SE = fitsread([path_data,cell2mat(dataID(k))]);
% long exposure images
im_LE = sum(im_SE,3);

%% TEST DETECTION


