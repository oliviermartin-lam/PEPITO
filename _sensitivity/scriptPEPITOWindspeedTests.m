clear all;close all;
savePath = '/home/omartin/Projects/PEPITO/Simulations/';
%% LOAD THE DATA
if false
    im = fitsread([savePath,'fits/simulatedPSFs/simulatedLEimagesVwindspeed.fits']);
   
    
    im2 = fitsread([savePath,'fits/simulatedPSFs/simulatedSEimagesVwindspeed.fits']);
    
    im0 = zeros(npsf,npsf,nSrc,nW,6);
    im0(:,:,:,:,1:5) = im2;
    im0(:,:,:,:,6) = im;
    
    wVect = [5,10,15,20,25,30,40,50];
    
else
    im0 = fitsread([savePath,'fits/simulatedPSFs/simulatedSEimagesVwindspeed50ms.fits']);      
    wVect = [50];
end

npsf = size(im0,1);
nIm = size(im0,3);
nW   = size(im0,4);
%% Define a telescope
D      = 1;
fov    = 0.2;
resCam = 60;
resPup = 60;
Fe     = 1000;
cobs   = 0.39;
tel    = telescope(D,'fieldOfViewInArcMin',fov,'resolution',resPup,...
    'samplingTime',1/Fe,'obstructionRatio',cobs);
% Define sources
n    = 2;
zS   = linspace(-fov*60/2,fov*60/2,n);             
aS   = pi/4*ones(1,n);
ngs  = source('wavelength',photometry.V0,'zenith',zS*constants.arcsec2radian,'azimuth',aS);
nSrc = numel(ngs);
atm.wavelength = ngs(1).photometry;
loD     = constants.radian2mas*ngs(1).wavelength/tel.D;
psInMas = loD/2;
%% DEFINE ATMOSPHERE
iW = 1;

wS = wVect(iW);
L0 = 25;
alt= 16e3;
w0 = 1;
wD = pi/4;
r0 = 0.1;
wS = 15;
atm = atmosphere(photometry.V0,r0,L0,'altitude',alt,'fractionnalR0',w0,...
    'windSpeed',wS,'windDirection',wD);
% Define otf0
psf0= im0(:,:,1,iW,1);
otf0= fourierTools.psf2otf(psf0);
otf0=otf0/max(otf0(:));
psf2= im0(:,:,2,iW,1);
otf2= fourierTools.psf2otf(psf2);
otf2=otf2/max(otf2(:));
ATFsimu = otf2./otf0;

%COMPUTING ATF
dk    = 11;
validPhase = true(dk);
otfDL = fourierTools.pupil2otf(tel.pupil,0*tel.pupil,1);
otfDL = psfrTools.interpolateOtf(otfDL,size(otf0,1));
otfDL = otfDL/max(otfDL(:));
mskOtf= otfDL>1e-8;


%% ----------------------- PSF MODEL VERIFICATION ---------------------- %
close all;
fvu = zeros(2,nIm-1);

% Get the diffraction-limit OTF
oDL = psfrTools.zonalCovarianceToOtf(zeros(dk^2),size(otf0,1),tel.D,tel.D/(dk-1),validPhase);
msk = oDL/max(oDL(:)) > 1e-5;
    
% ---------- Get the simulated ATF
%ATFsimu = mskOtf.*real(fourierTools.psf2otf(im(:,:,2,iW))./fourierTools.psf2otf(im(:,:,1,iW)));
%ATFsimu(ATFsimu<0)=0;

% ---------- Get the analytical ATF
% Covariance matrices
[Css,Cgs]= phaseStats.spatioAngularCovarianceMatrix(dk,tel.D,atm,ngs(2),'srcCC',ngs(1));
% TT filter
[X,Y]    = meshgrid((-1+1/dk:2/dk:1-1/dk),(-1+1/dk:2/dk:1-1/dk));
TT       = [X(:),Y(:)];
FTT      = TT*pinv(TT);
% TT anisoplanatism cov matrix
CaniTT   = FTT*(2*Css - Cgs{1}' - Cgs{1}')*FTT'/sqrt(2);
tmp      = psfrTools.zonalCovarianceToOtf(CaniTT,size(otf0,1),tel.D,tel.D/(dk-1),validPhase);
tmp(msk) = tmp(msk)./oDL(msk);
% Analytical ATF
ATFTT  = mskOtf.*tmp/tmp(end/2+1,end/2+1);
fvu(1) = 1e2*psfrTools.getFVU(ATFsimu,ATFTT);

% ---------- Get the Analytical PSF
psfMod = fourierTools.otf2psf(otf0.*ATFTT);
psfMod = psfMod/sum(psfMod(:));
imi    = abs(im0(:,:,2,iW,1)/sum(sum(im0(:,:,2,iW,1))));
fvu(2) = 1e2*psfrTools.getFVU(imi,psfMod);
        
% Display
idx   = npsf:npsf-1:npsf*(npsf-1)+1;
aSimu = abs(ATFsimu(:));
aTT   = ATFTT(:);
u     = linspace(-1,1,npsf);
pSimu = imi(:);
pTT   = psfMod(:);
x     = psInMas*1e-3*linspace(-npsf/2,npsf/2,npsf);%linspace(-npsf/2+1/npsf,npsf/2-1/npsf);

    
% PSF cutting plot
h=figure;
subplot(1,3,1)
plot(x,log10(diag(imi)),'b-');hold on;
plot(x,log10(diag(psfMod)),'r-.');
plot(x,log10(abs(diag(imi-psfMod))),'k:');
xlabel('Off-axis distance (arcsec)','interpreter','latex','FontSize',24);
ylabel('PSF intensity','interpreter','latex','FontSize',24);
title('Main diagonal','interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
ylim([-8,-1.5]);
subplot(1,3,2)
plot(x,log10(pSimu(idx)),'b-');hold on;
plot(x,log10(pTT(idx)),'r-.');
plot(x,log10(abs(pTT(idx)-pSimu(idx))),'k:');
xlabel('Off-axis distance (arcsec)','interpreter','latex','FontSize',24);
title('Secondary diagonal','interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
ylim([-8,-1.5]);
subplot(1,3,3)
plot(x(end/2+1:end),log10(radial(imi)),'b-');hold on;
plot(x(end/2+1:end),log10(radial(psfMod)),'r-.');
plot(x(end/2+1:end),log10(radial(abs(imi-psfMod))),'k:');
xlabel('Off-axis distance (arcsec)','interpreter','latex','FontSize',24);
title('Azimuthal average','interpreter','latex','FontSize',24);
legend({'Simulation','Model','Residual'},'interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
ylim([-8,-1.5]);
    
% PSF 2D maps
h=figure;
subplot(1,3,1);
imagesc(log10(imi),[-5,-1.2]);
title('Simulated PSF','interpreter','latex','FontSize',24)
set(gca,'Xtick',[],'Ytick',[]);
subplot(1,3,2);
imagesc(log10(psfMod),[-5,-1.2]);
set(gca,'Xtick',[],'Ytick',[]);
title('PSF Model','interpreter','latex','FontSize',24)
subplot(1,3,3);
imagesc(log10(abs(imi-psfMod)),[-5,-1.2]);
set(gca,'Xtick',[],'Ytick',[]);
title('Residual','interpreter','latex','FontSize',24)
    

%% ------------------------ CN2 RETRIEVAL ------------------------------ %%



% Instantiate the PRIMEclass
Wn = zeros(nW,1);
fW = zeros(nW,1);
Ln = zeros(nW,1);
fL = zeros(nW,1);


for iW=1:nW
    
    otfC= fourierTools.psf2otf(im0(:,:,1,iW));
    otfC=otfC/max(otfC(:));
    %imLE = squeeze(mean(im0(:,:,2,iW,:),5));

    imLE = im0(:,:,2);
    
    prime   = primeClass(tel,ngs(2),ngs(1),atm,imLE/sum(imLE(:)),psInMas,otfC,...
    [],[],'nLayer',1,'weightInit',atm.r0^(-5/3),'fitWeights',true,'heightInit',16e3,...
    'fitHeight',false,'nptZonal',dk,'idxPhase',validPhase,...
    'flagAnisokinetism',true,'fitOuterScale',1);
    % Cn2h retrieval
    prime  = prime.getParameters();
    Wn(iW) = prime.xfinal(1);
    fW(iW) = prime.dweights;
    Ln(iW) = prime.xfinal(2);
    fL(iW) = prime.dL0;
end

wRef = 45.5947;%atm.r0^(-5/3);
dWn = 100*bsxfun(@rdivide,bsxfun(@minus,Wn,wRef'),wRef');
%fitswrite([dWn,fW],[savePath,'fits/cn2hErrorVWindspeed.fits']);



h = figure;
errorbar(wVect,(dWn),fW*3,'bs--','MarkerFaceColor','b','MarkerSize',7);
hold on;
xlabel('Windspeed value (m/s)','interpreter','latex','FontSize',24);
ylabel('$C_n^2(h)$ error (\%)','interpreter','latex','FontSize',24);
psfrTools.makeAxisSquared(h);
xlim([0,60]);
plot([12.5,12.5],[ylim()],'k--');
plot([35,35],[ylim()],'k--');
plot([xlim()],[0,0],'k-.');
text(1,23,'Convergence','interpreter','latex','FontSize',18);
text(17,23,'Optimal region','interpreter','latex','FontSize',18);
text(40,23,'Finite exposure','interpreter','latex','FontSize',18);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');




