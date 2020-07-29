clear all;close all;
savePath = '/home/omartin/Projects/PEPITO/Simulations/';
%% LOAD THE DATA
atmCase = 1;
    
if atmCase == 1
    im = fitsread([savePath,'fits/simulatedPSFs/simulatedLongExposureImages_recentered_dayTime.fits']);
    r0 = 0.06;
else
    im = fitsread([savePath,'fits/simulatedPSFs/simulatedLongExposureImages_recentered_dayTime_case2.fits']);
    r0 = 0.1;
end
npsf = size(im,1);
nSrc = size(im,3);


%Get FWHM/Aspect ratio
FWHM = zeros(1,nSrc);
aspR = zeros(1,nSrc);
for iSrc=1:nSrc
    [F1,F2,e] = psfrTools.getFWHM(im(:,:,iSrc),1,4);
    FWHM(iSrc) = max(F1,F2);
    aspR(iSrc) = e;
end

    
%% DEFINE ATMOSPHERE
L0 = 25;
alt= [1,4,8,16]*1e3;
w0 = [50,25,10,15]/100;
wS = [5,25,15,10];
wD = [0,pi/4,pi/2,pi];
atm = atmosphere(photometry.V0,r0,L0,'altitude',alt,'fractionnalR0',w0,...
    'windSpeed',wS,'windDirection',wD);
%% Define the telescope/sources
D      = 1.3;
fov    = 5;
resCam = 30;
nRes   = 30;
Fe     = 100;
cobs   = 0.39;
tel    = telescope(D,'fieldOfViewInArcMin',fov,'resolution',nRes,...
    'samplingTime',1/Fe,'obstructionRatio',cobs);

% Define sources
n       = 31;
zS      = linspace(-fov*60/2,fov*60/2,n);               %
aS      = pi/4*ones(1,n);
ngs     = source('wavelength',photometry.Kp2,'zenith',zS*constants.arcsec2radian,'azimuth',aS);
loD     = constants.radian2mas*ngs(1).wavelength/tel.D;
psInMas = loD/2;

atm.wavelength = ngs(1).photometry;

% Define otf0
psf0= im(:,:,1);
otf0= fourierTools.psf2otf(psf0);
otf0=otf0/max(otf0(:));

%COMPUTING ATF
dk    = 11;
validPhase = true(dk);
otfDL = fourierTools.pupil2otf(tel.pupil,0*tel.pupil,1);
otfDL = psfrTools.interpolateOtf(otfDL,size(otf0,1));
otfDL = otfDL/max(otfDL(:));
mskOtf= otfDL>1e-8;


%% ----------------------- PSF MODEL VERIFICATION ---------------------- %
close all;
nIm = size(im,3);
fvu = zeros(2,nIm-1);

% Get the diffraction-limit OTF
oDL = psfrTools.zonalCovarianceToOtf(zeros(dk^2),size(otf0,1),tel.D,tel.D/(dk-1),validPhase);
msk = oDL/max(oDL(:)) > 1e-5;
    
for iSrc=1:nIm
    % ---------- Get the simulated ATF
    ATFsimu = mskOtf.*real(fourierTools.psf2otf(im(:,:,iSrc))./fourierTools.psf2otf(im(:,:,1)));
    ATFsimu(ATFsimu<0)=0;
    
    % ---------- Get the analytical ATF 
    % Covariance matrices
    [Css,Cgs]   = phaseStats.spatioAngularCovarianceMatrix(dk,tel.D,atm,ngs(iSrc),'srcCC',ngs(1));    
    % TT filter
    [X,Y]       = meshgrid((-1+1/dk:2/dk:1-1/dk),(-1+1/dk:2/dk:1-1/dk));
    TT          = [X(:),Y(:)];
    FAA         = TT*pinv(TT);
    % TT anisoplanatism cov matrix
    CaniTT      = FAA*(2*Css - Cgs{1}' - Cgs{1}')*FAA'/sqrt(2);
    tmp         = psfrTools.zonalCovarianceToOtf(CaniTT,size(otf0,1),tel.D,tel.D/(dk-1),validPhase);
    tmp(msk)    = tmp(msk)./oDL(msk);
    % Analytical ATF
    ATFTT       = mskOtf.*tmp/tmp(end/2+1,end/2+1);    
    fvu(1,iSrc-1) = 1e2*psfrTools.getFVU(ATFsimu,ATFTT);
    
    % ---------- Get the Analytical PSF
    psfMod        = fourierTools.otf2psf(otf0.*ATFTT);
    psfMod        = psfMod/sum(psfMod(:));
    imi           = abs(im(:,:,iSrc)/sum(sum(im(:,:,iSrc))));
    fvu(2,iSrc-1) = 1e2*psfrTools.getFVU(imi,psfMod);
    
    
    if iSrc ==16
        % Display
        idx   = npsf:npsf-1:npsf*(npsf-1)+1;
        aSimu = abs(ATFsimu(:));
        aTT   = ATFTT(:);
        u     = linspace(-1,1,npsf);
        pSimu = imi(:);
        pTT   = psfMod(:);
        x     = psInMas*1e-3*linspace(-npsf/2,npsf/2,npsf);%linspace(-npsf/2+1/npsf,npsf/2-1/npsf);
        
        % ATF cutting plots
        h = figure;
        subplot(1,3,1)
        plot(u,aSimu(idx),'b-');hold on;
        plot(u,aTT(idx),'r-.');
        plot(u,abs(aTT(idx)-aSimu(idx)),'k:');
        xlabel('$D/\lambda$ units','interpreter','latex','FontSize',24);
        ylabel('Spatial frequency filter','interpreter','latex','FontSize',24);
        title('Main diagonal','interpreter','latex','FontSize',24);
        set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
        ylim([0,1]);
        
        subplot(1,3,2)
        plot(u,diag(ATFsimu),'b-');hold on;
        plot(u,diag(ATFTT),'r-.');
        plot(u,abs(diag(ATFTT-ATFsimu)),'k:');
        xlabel('$D/\lambda$ units','interpreter','latex','FontSize',24);
        title('Secondary diagonal','interpreter','latex','FontSize',24);
        set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
        ylim([0,1]);
        
        subplot(1,3,3)
        plot(u(end/2+1:end),radial(ATFsimu),'b-');hold on;
        plot(u(end/2+1:end),radial(ATFTT),'r-.');
        plot(u(end/2+1:end),radial(abs(ATFTT-ATFsimu)),'k:');
        xlabel('$D/\lambda$ units','interpreter','latex','FontSize',24);
        title('Azimuthal average','interpreter','latex','FontSize',24);
        legend({'Simulation','Model','Residual'},'interpreter','latex','FontSize',24);
        set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
        ylim([0,1]);
        
        % ATF 2D maps
        h=figure;
        subplot(1,3,1);
        imagesc(ATFsimu,[0,1]);
        title('Simulation','interpreter','latex','FontSize',24)
        set(gca,'Xtick',[],'Ytick',[]);
        subplot(1,3,2);
        imagesc(ATFTT,[0,1]);
        set(gca,'Xtick',[],'Ytick',[]);
        title('Model','interpreter','latex','FontSize',24)
        subplot(1,3,3);
        imagesc(abs(ATFsimu-ATFTT),[0,1]);
        set(gca,'Xtick',[],'Ytick',[]);
        title('Residual','interpreter','latex','FontSize',24)
        
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
    end
end

h = figure;
plot(zS(2:end),fvu(2,:),'bo--','MarkerFaceColor','b','MarkerSize',7);
ylabel('FVU on PSF (\%)','interpreter','latex','FontSize',24);
xlabel('Off-axis distance (arcsec)','interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
psfrTools.makeAxisSquared(h);


%% ------------------------ CN2 RETRIEVAL ------------------------------ %%
% CCD config
tExp      = 1e-2;
nExp      = 120/tExp;
ron       = 40;
darkBg    = 0.05;
eCCD      = sqrt(nExp*(ron^2 + (darkBg*tExp)^2));
QE        = 0.7;
throughput= 0.9;
deltaL    = 0.1*350; %10nm-> Angstrom
F0        = 43.6*1e6*deltaL;
F2DN      = tExp*nExp*throughput*QE*tel.area;

% Instantiation
mag       = 0:8;
nMag      = length(mag);
nIter     = 10;
Wn        = zeros(4,nMag,nIter);
fW        = 0*Wn;

% Initial guess
h0 = [atm.layer.altitude];
w0 = [2,1.5,.5,1];

% Binning the image
idN   = [2,5,10,15,20,25];
nRef  = length(idN);
ps    = constants.radian2mas*ngs(1).wavelength/tel.D/2;
nBin  = round(60*ps/500);
psBin = npsf*ps/nBin;
skyBg = (psBin*1e-3)^2*10^(-0.4*6.5)*F0*F2DN;

imN = zeros(nBin,nBin,nRef);
for iSrc=1:nRef
    imN(:,:,iSrc) = psfrTools.interpolateOtf(im(:,:,idN(iSrc)),nBin);
end
imRef = reshape(imN,nBin,nBin*nRef);



% Instantiate the PRIMEclass
prime   = primeClass(tel,ngs(idN),ngs(1),atm,imRef,psBin,otf0,...
    [],[],'nLayer',4,'weightInit',w0,'fitWeights',true,'heightInit',h0,...
    'fitHeight',false,'nptZonal',dk,'idxPhase',validPhase,...
    'flagAnisokinetism',true);

for k=1:nMag
    for t=1:nIter
        % Noisify the image
        for iSrc = 1:nRef
            idx = 1  +(iSrc-1)*nBin:iSrc*nBin;
            tmp = imRef(:,idx)*10^(-0.4*mag(k))*F0*F2DN/sum(sum(imRef(:,idx)));
            prime.psfSky(:,idx) = poissrnd(tmp + skyBg) - skyBg + eCCD*randn(nBin);
        end
        % Cn2h retrieval
        prime       = prime.getParameters();
        Wn(:,k,t) = prime.xfinal;
        fW(:,k,t) = prime.dweights;
    end
end

wRef = atm.r0^(-5/3)*[atm.layer.fractionnalR0];
wRef = [3.1912,1.5095,0.8584,0.7996];
dWn = 100*bsxfun(@rdivide,bsxfun(@minus,Wn,wRef'),wRef');    
fitswrite([dWn,fW],[savePath,'fits/cn2hErrorVMag_DayTime.fits']);

    
% Plot
dW = mean(abs(dWn),3);
eW = sqrt(var(dWn,[],3) + mean(fW.^2,3))/sqrt(nIter);
dW1= squeeze(dW(1,:));
eW1= squeeze(eW(1,:));
dW2= squeeze(dW(2,:));
eW2= squeeze(eW(2,:));
dW3= squeeze(dW(3,:));
eW3= squeeze(eW(3,:));
dW4= squeeze(dW(4,:));
eW4= squeeze(eW(4,:));


h = figure;
errorbar(mag,dW1,eW1,'ks--','MarkerFaceColor','k','MarkerSize',7,'LineWidth',1.2);hold on;
errorbar(mag,dW2,eW2,'bo--','MarkerFaceColor','b','MarkerSize',7,'LineWidth',1.2);
errorbar(mag,dW3,eW3,'rd--','MarkerFaceColor','r','MarkerSize',7,'LineWidth',1.2);
errorbar(mag,dW4,eW4,'mp--','MarkerFaceColor','m','MarkerSize',7,'LineWidth',1.2);
xlabel('Calibration star magnitude (mag)','interpreter','latex','FontSize',24);
ylabel('$C_n^2(h)$ error (\%)','interpreter','latex','FontSize',24);
legend({'1km','4km','8km','16km'},...
    'interpreter','latex','FontSize',24,'Location','northwest','AutoUpdate','off');
psfrTools.makeAxisSquared(h);
xlim([mag(1)-1,mag(end)+1]);
plot([xlim()],[1,1],'k--');
text(7.5,1.2,'1\%','interpreter','latex','FontSize',24);
plot([xlim()],[10,10],'k--');
text(7.5,11.8,'10\%','interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
set(gca,'YScale','log')

%r0ref = 0.3543
%0.5,0.25,0.1,0.15

