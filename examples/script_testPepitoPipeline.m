clear all;
close all;
% Get simulated PSFs
a  = fitsread('simulatedShortExposureImages.fits');
% Image format
nBox   = size(a,1);
nStars = size(a,3);
% Photometric characteristics
wvl      = 500e-9;
D        = 0.5;
obs      = 0.3;
magStars = [7.13 8.3 11.25 10.29 10.36 10.37 12.14];
tExp     = 10e-3;
nExp     = 3000;%size(a,4);
zeroPoint= 995.5*1e6;%photon/m2/s/Angstrom
lWidth   = 50*0.1;
QE       = .8;
throughput= .7;
Flux2ADU = zeroPoint*D^2*(1-obs)^2*tExp*lWidth*QE*throughput;
fluxS    = 10.^(-0.4*magStars)*Flux2ADU;
% Noise
psInArcsec=180*3600/pi*wvl/D/2;
ron      = 1.;
darkBg   = 0.1;
eN       = sqrt((ron^2 + (darkBg*tExp)^2));
skyBg    = (psInArcsec*1e-3)^2*10^(-0.4*24)*Flux2ADU;

redo = true;
flagNoise = false;


xSinPix  = [45.3 57.1  98.2  9.05 24.32 3.09 78.2]+15.5;
ySinPix  = [2.3  72.56 56.21 42.2 90.1  5.43 27.2]+15.5;
subX     = (xSinPix - floor(xSinPix));
subY     = (ySinPix - floor(ySinPix));
xSinPix  = round(xSinPix);
ySinPix  = round(ySinPix);

    
if redo    
    %Preparing the loop ...
    maxRes  = max(max(xSinPix),max(ySinPix)) + nBox/2;
    idxi    = xSinPix - nBox/2 + 1;
    idxf    = xSinPix + nBox/2;
    idyi    = ySinPix - nBox/2 + 1;
    idyf    = ySinPix + nBox/2;
    
    % Loop on stars
    im_    = zeros(maxRes,maxRes,nExp);
    for k=1:nExp
        for iS = 1:nStars
            % Flux scaling
            psf_i = a(:,:,iS,k);
            psf_i = psf_i/sum(psf_i(:));
            % Translate the image
            psf_i = pepito.translateImage(psf_i,[subX(iS),subY(iS)]);
            psf_i(psf_i<0) = 0;
            % Update the image
            idX          = round(idxi(iS):idxf(iS));
            idY          = round(idyi(iS):idyf(iS));
            im_(idX,idY,k) = im_(idX,idY,k) + psf_i*fluxS(iS);
        end
	if flagNoise
        	im_ (:,:,k) = poissrnd(im_(:,:,k) + skyBg)  -skyBg + eN*randn(maxRes);
	end
    end
    
    fitswrite(im_,'fakeImage_noNoise.fits');
else       
    im_ = fitsread('fakeImage_noNoise.fits');
end

close all;
figure;
imagesc((abs(sum(im_,3))))
            

%% ---------------------------- RUN PEPITO ----------------------------- %%
psInMas= psInArcsec*1e3;
RA     = [0,11.17,135.15,172.07,191.7,247.15,248.32];
DEC    = zeros(1,numel(RA));
nBox   = 30;
nLayer = 7;
hInit  = [0,.5,1,2,4,8,16]*1e3;
wInit  = [40,20,5,7,15,5,8]/100*0.1^(-5/3);

pep   = pepito(im_,psInMas,D,RA,DEC,nBox,nLayer,'npt',21,'wvl',500e-9,'nref',1, ...
    'weightInit',wInit,'heightInit',hInit);


pep = pep.runPepito([xSinPix,ySinPix,fluxS]);

  
