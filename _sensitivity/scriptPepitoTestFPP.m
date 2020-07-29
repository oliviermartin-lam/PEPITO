clear all;close all;
savePath = '/home/omartin/Projects/PEPITO/Simulations/';
%% LOAD THE DATA
im = fitsread([savePath,'fits/simulatedPSFs/simulatedLongExposureImages_recentered_multilayers.fits']);
npsf = size(im,1);
nSrc = size(im,3);

%% --------------- Section 3.1: Digital anisokinetism ------------------ %%

tmp = fitsread([savePath,'fits/LEPSFstat_50cm.fits']);  
F50  = tmp(1:end/2);
e50  = tmp(end/2+1:end);
tmp = fitsread([savePath,'fits/LEPSFstat_1m.fits']);  
F1  = tmp(1:end/2);
e1  = tmp(end/2+1:end);
zS  = linspace(0,30,length(F50));
h = figure;
plot(zS,F50,'bo--','MarkerFaceColor','b','MarkerSize',7);
hold on;
plot(zS,F1,'rs-.','MarkerFaceColor','r','MarkerSize',7);
ylabel('PSF FWHM (mas)','interpreter','latex','FontSize',24);
xlabel('Off-axis distance (arcsec)','interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
legend({'50cm-telescope','1m-telescope'},'interpreter','latex','FontSize',24,'Location','northwest');
psfrTools.makeAxisSquared(h);

h = figure;
plot(zS,e50,'bo--','MarkerFaceColor','b','MarkerSize',7);
hold on;
plot(zS,e1,'rs-.','MarkerFaceColor','r','MarkerSize',7);
ylabel('PSF Aspect ratio','interpreter','latex','FontSize',24);
xlabel('Off-axis distance (arcsec)','interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
legend({'50cm-telescope','1m-telescope'},'interpreter','latex','FontSize',24,'Location','northwest');
psfrTools.makeAxisSquared(h);

%% ------------------ Section 3.2: Model validation -------------------- %%

singleLayer = false;
% single layer case
if singleLayer
    r0 = 0.1;
    L0 = 25;
    w0 = 0.5*atm.r0^(-5/3);
    h0 = 8e3;
    nL = 1;
    atm = atmosphere(photometry.V0,r0,L0,'altitude',h0,'fractionnalR0',1,...
    'windSpeed',0,'windDirection',0);
else
    % multi-layers case
    L0 = 25;
    r0 = 0.15*cos(30*pi/180)^(3/5);
    w0 = [0.6 0.3 0.1];
    h0 = [4 8 16]*1e3;
    nL = 3;
    atm = atmosphere(photometry.V0,r0,L0,'altitude',h0,'fractionnalR0',w0,...
    'windSpeed',[5.5,9.5,6.3],'windDirection',[-19.57,-60.82,-10.4]*pi/180);
end

% Define a telescope
D   = 0.5;
cobs= 0.39;
%D   = 1;
%cobs = 0.3;

fov = 1;
nRes= 100;
Fe  = 1000;

tel = telescope(D,'fieldOfViewInArcMin',fov,'resolution',nRes,...
    'samplingTime',1/Fe,'obstructionRatio',cobs);

% Define sources
n         = 16;
zS       = linspace(0,fov*60/2,n)*constants.arcsec2radian;
aS       = pi/4*ones(1,n);
ngs     = source('wavelength',photometry.V0,'zenith',zS,'azimuth',aS);
nSrc    = numel(ngs);
psInMas = constants.radian2mas*ngs(1).wavelength/tel.D/2;

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


close all;
nIm = size(im,3);
fvu = zeros(2,nIm-1);

for iSrc=2:nIm
    ATFsimu = mskOtf.*real(fourierTools.psf2otf(im(:,:,iSrc))./fourierTools.psf2otf(im(:,:,1)));
    ATFsimu(ATFsimu<0)=0;
    % Model 1: AA filtering on spatial phase covariance
    
    oDL         = psfrTools.zonalCovarianceToOtf(zeros(dk^2),size(otf0,1),tel.D,tel.D/(dk-1),validPhase);
    msk         = oDL/max(oDL(:)) > 1e-5;
    [Css,Cgs]   = phaseStats.spatioAngularCovarianceMatrix(dk,tel.D,atm,ngs(iSrc),'srcCC',ngs(1));
    
    % AA filter
    [X,Y]       = meshgrid((-1+1/dk:2/dk:1-1/dk),(-1+1/dk:2/dk:1-1/dk));
    AA          = [X(:),Y(:)];
    FAA         = AA*pinv(AA);
    % AA cov matrix
    CaniAA      = FAA*(2*Css - Cgs{1}' - Cgs{1}')*FAA'/sqrt(2);%1.1855^2;
    tmp         = psfrTools.zonalCovarianceToOtf(CaniAA,size(otf0,1),tel.D,tel.D/(dk-1),validPhase);
    tmp(msk)    = tmp(msk)./oDL(msk);
    % ATF
    ATFAA       = mskOtf.*tmp/tmp(end/2+1,end/2+1);
    
    fvu(1,iSrc-1) = 1e2*psfrTools.getFVU(ATFsimu,ATFAA);
    psfMod        = fourierTools.otf2psf(otf0.*ATFAA);
    psfMod        = psfMod/sum(psfMod(:));
    imi           = abs(im(:,:,iSrc)/sum(sum(im(:,:,iSrc))));
    fvu(2,iSrc-1) = 1e2*psfrTools.getFVU(imi,psfMod);

    if iSrc ==16
        % Display
        idx   = npsf:npsf-1:npsf*(npsf-1)+1;
        aSimu = abs(ATFsimu(:));
        aAA   = ATFAA(:);
        u     = linspace(-1,1,npsf);
        pSimu = imi(:);
        pAA   = psfMod(:);
        x     = psInMas*1e-3*linspace(-npsf/2+1/npsf,npsf/2-1/npsf);
        
        % ATF cutting plots
        h = figure;
        subplot(1,3,1)
        plot(u,aSimu(idx),'b-');hold on;
        plot(u,aAA(idx),'r-.');
        plot(u,abs(aAA(idx)-aSimu(idx)),'k:');
        xlabel('$D/\lambda$ units','interpreter','latex','FontSize',24);
        ylabel('Spatial frequency filter','interpreter','latex','FontSize',24);
        title('Main diagonal','interpreter','latex','FontSize',24);
        set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
        ylim([0,1]);
        
        subplot(1,3,2)
        plot(u,diag(ATFsimu),'b-');hold on;
        plot(u,diag(ATFAA),'r-.');
        plot(u,abs(diag(ATFAA-ATFsimu)),'k:');
        xlabel('$D/\lambda$ units','interpreter','latex','FontSize',24);
        title('Secondary diagonal','interpreter','latex','FontSize',24);
        set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
        ylim([0,1]);
        
        subplot(1,3,3)
        plot(u(end/2+1:end),radial(ATFsimu),'b-');hold on;
        plot(u(end/2+1:end),radial(ATFAA),'r-.');
        plot(u(end/2+1:end),radial(abs(ATFAA-ATFsimu)),'k:');
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
        imagesc(ATFAA,[0,1]);
        set(gca,'Xtick',[],'Ytick',[]);
        title('Model','interpreter','latex','FontSize',24)
        subplot(1,3,3);
        imagesc(abs(ATFsimu-ATFAA),[0,1]);
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
        plot(x,log10(pAA(idx)),'r-.');
        plot(x,log10(abs(pAA(idx)-pSimu(idx))),'k:');
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
plot(zS(2:end)*constants.radian2arcsec,fvu(1,:),'bo','MarkerFaceColor','b','MarkerSize',7);
ylabel('FVU on ATF (\%)','interpreter','latex','FontSize',24);
xlabel('Off-axis distance (arcsec)','interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
psfrTools.makeAxisSquared(h);


fvu1m = fitsread([savePath,'fits/fvuPSF_1mVis.fits']);
h = figure;
plot(zS(2:end)*constants.radian2arcsec,fvu(2,:),'bo--','MarkerFaceColor','b','MarkerSize',7);
hold on;
plot(zS(2:end)*constants.radian2arcsec,fvu1m,'rs-.','MarkerFaceColor','r','MarkerSize',7);
ylabel('FVU on PSF (\%)','interpreter','latex','FontSize',24);
xlabel('Off-axis distance (arcsec)','interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
legend({'50cm-telescope','1m-telescope'},'interpreter','latex','FontSize',24,'Location','northwest');
psfrTools.makeAxisSquared(h);


%% ---------------------- Section 3.3: Performance ---------------------- %%
zenithAngle = 45;
airmass     = 1/cos(zenithAngle*pi/180);
r0          = 0.157*airmass^(-3/5);
atm2 = atmosphere(photometry.V0,r0,L0,'altitude',8e3,'fractionnalR0',1,...
    'windSpeed',0,'windDirection',0);

src2 = source('zenith',12*constants.arcsec2radian,'azimuth',pi/4);
redo = false;

if redo
    nT = 500;
    WL = zeros(2,nT);
    fW = 0*WL;
    t0 = zeros(1,nT);
    
    fppr = focalPlaneProfiler(tel,src2,ngs(1),atm2,otf0,psInMas,otf0,...
        [],[],'nLayer',1,'weightInit',r0^(-5/3),'fitWeights',true,'heightInit',8e3,...
        'fitHeight',false,'fitOuterScale',1,'nptZonal',dk,'idxPhase',validPhase,'flagAnisokinetism',true);
    psfModel       = fppr.psfModel([r0^(-5/3),L0],fppr.otfKernel);
    fppr.psfSky     = psfModel;
    s = 5;
    
    for t=1:nT
        tic;
        fppr.L0Init    = L0+trandn((20-L0)/s,(30-L0)/s)*s;
        fppr.weightInit= r0^(-5/3);
        fppr           = fppr.getParameters();
        WL(:,t)       = fppr.xfinal;
        fW(:,t)       = [fppr.dweights fppr.dL0];
        t0(t) = toc();
    end
    fitswrite([WL,fW],[savePath,'fits/Cn2h_L0_ErrorVinitGuess.fits']);
    
    WC  = WL(1:end/2,:);
    fC  = fW(1:end/2,:);
    WL  = WL(end/2+1:end,:);
    fL  = fW(end/2+1:end,:);
else
    out = fitsread([savePath,'fits/Cn2h_L0_ErrorVinitGuess.fits']);  
    WC  = out(1:end/2,1:end/2);
    fC  = out(1:end/2,end/2+1:end);
    WL  = out(end/2+1:end,1:end/2);
    fL  = out(end/2+1:end,end/2+1:end);
    nT  = size(WL);
end


dL = 1e2*(WL-L0)/L0;
dC = 1e2*(WC-r0^(-5/3))/r0^(-5/3);

id = abs(dL)<5;
h = figure;
histogram(dC(id),30,'Normalization','probability');hold on;
histogram(dL(id),30,'Normalization','probability');
ylabel('Probability','interpreter','latex','FontSize',24);
xlabel('Estimation error (\%)','interpreter','latex','FontSize',24);
psfrTools.makeAxisSquared(h);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
legend({'$C_n^2(h)$','$L_0$'},'interpreter','latex',...
    'FontSize',24,'Location','northeast','AutoUpdate','off');

mean(dL(id)),std(dL(id))
mean(dC(id)),std(dC(id))


%% ---------------------- Section 4.1: Cn2h = f(t) ---------------------- %%
close all;

redo = false;

 nF = (100:100:5000);
 r0 = 0.15*cos(30*pi/180)^(3/5);
 w0 = r0^(-5/3);
    
 if redo
     D   = 1;
     %Load images
     if D==0.5
         cobs = 0.39;
         zS   = 12;
         im   = fitsread([savePath,'fits/simulatedPSFs/simulatedShortExposureImages_recentered_singlelayer_50cmVis.fits']);
     elseif D==1
         cobs = 0.3;
         im   = fitsread([savePath,'fits/simulatedPSFs/simulatedShortExposureImages_recentered_singlelayer_1mVis.fits']);
         zS   = 19;
     end
     
     % Telescope config
     tel  = telescope(D,'fieldOfViewInArcMin',1,'resolution',100,...
         'samplingTime',1e-3,'obstructionRatio',cobs);
     dk   = 11;
     validPhase = true(dk);
     % Atmosphere config
     atm2 = atmosphere(photometry.V0,r0,25,'altitude',8e3,'fractionnalR0',1,...
         'windSpeed',0,'windDirection',0);
     % Detector config
     npsf = size(im,1);
     ps   = constants.radian2mas*500e-9/D/2;
     src2 = source('wavelength',photometry.V0,'azimuth',pi/4,'zenith',zS*constants.arcsec2radian);
     % Instantiating outputs
     
     nT = length(nF);
     cn2h = zeros(1,nT);
     pcn2h = zeros(1,nT);
     %Loop on acquisition time
     for t=1:nT
         im0  = squeeze(mean(im(:,:,1,1:nF(t)),4));
         im2  = squeeze(mean(im(:,:,2,1:nF(t)),4));
         otf0 = fourierTools.psf2otf(im0);
         
         fppr = focalPlaneProfiler(tel,src2,source('wavelength',photometry.V0),atm2,im2,ps,otf0,...
             [],[],'nLayer',1,'weightInit',w0,'fitWeights',true,'heightInit',8e3,...
             'fitHeight',false,'nptZonal',dk,'idxPhase',validPhase,'flagAnisokinetism',true);
         
         fppr = fppr.getParameters();
         cn2h(t) = fppr.xfinal;
         pcn2h(t) = fppr.dweights;
     end
     
     fitswrite([cn2h,pcn2h],[savePath,'fits/cn2hErrorVtime_1m.fits']);
 else
     tmp    = fitsread([savePath,'fits/cn2hErrorVtime_50cm.fits']);
     cn2h50 = tmp(1:end/2);
     pcn2h50= tmp(end/2+1:end);
     tmp    = fitsread([savePath,'fits/cn2hErrorVtime_1m.fits']);
     cn2h1  = tmp(1:end/2);
     pcn2h1 = tmp(end/2+1:end);
 end
%
dt= 1e-2*nF(2:end);
ns = 20;
dW50 = 1e2*smooth(abs(cn2h50(2:end) - w0)/w0,ns)';
dW1  = 1e2*smooth(abs(cn2h1(2:end) - w0)/w0,ns)';
eW50 = 1e2*smooth(pcn2h50(2:end)/w0,ns)';
eW1  = 1e2*smooth(pcn2h1(2:end)/w0,ns)';

h = figure;
plot(dt,dW50,'bs-.','MarkerFaceColor','b','MarkerSize',5,'LineWidth',1.2);
hold on;
plot(dt,dW1,'ro--','MarkerFaceColor','r','MarkerSize',5,'LineWidth',1.2);
ylabel('$C_n^2(h)$ error (\%)','interpreter','latex','FontSize',24);
xlabel('Acquisition time (s)','interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
legend({'50cm-telescope','1m-telescope'},'interpreter','latex',...
    'FontSize',24,'Location','northeast','AutoUpdate','off');
psfrTools.makeAxisSquared(h);
inBetween = [(dW50-eW50) fliplr(dW50+eW50)];
fill([dt fliplr(dt)]', inBetween, 'b','FaceAlpha',.15,'LineStyle','none');
inBetween = [(dW1-eW1) fliplr(dW1+eW1)];
fill([dt fliplr(dt)]', inBetween, 'r','FaceAlpha',.15,'LineStyle','none');
plot([xlim()],[1,1],'k:','LineWidth',1.2)
text(5,2,'1\%-error','interpreter','latex','FontSize',24);


% tmp  = fitsread([savePath,'fits/FVU_aspectRatioVtime_50cm.fits']);  
% fvu50=1e2*sqrt(tmp(:,1:end/2));
% e50  =tmp(:,end/2+1:end);
% tmp  = fitsread([savePath,'fits/FVU_aspectRatioVtime_1m.fits']);  
% fvu1 =1e2*sqrt(tmp(:,1:end/2));
% e1   =tmp(:,end/2+1:end);
% 
% dt   = 1e-2*(500:250:5000);
% id50 = 1:16;
% id1  = 1:16;
% 
% 
% h = figure;
% plot(dt,mean(bsxfun(@rdivide,e50(id50,:),e50(id50,end)),1),'b-');
% hold on;
% plot(dt,mean(bsxfun(@rdivide,e1(id1,:),e1(id1,end)),1),'r--');
% ylabel('Normalized PSF aspect ratio','interpreter','latex','FontSize',24);
% xlabel('Acquisition time (s)','interpreter','latex','FontSize',24);
% set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
% legend({'50cm-telescope','1m-telescope'},'interpreter','latex',...
%     'FontSize',24,'Location','northeast','AutoUpdate','off');
% psfrTools.makeAxisSquared(h);
% plot([xlim()],[1,1],'k:','LineWidth',1.2)
% plot([xlim()],[1.01,1.01],'k-.','LineWidth',1.2)
% plot([xlim()],[.99,.99],'k--','LineWidth',1.2)
% text(5,1.003,'1\%-variation','interpreter','latex','FontSize',24);
% 
% h = figure;
% plot(dt,mean(fvu50(id50,:),1),'b-');
% hold on;
% plot(dt,mean(fvu1(id1,:),1),'r--');
% ylabel('FVU (\%)','interpreter','latex','FontSize',24);
% xlabel('Acquisition time (s)','interpreter','latex','FontSize',24);
% set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
% legend({'50cm-telescope','1m-telescope'},'interpreter','latex',...
%     'FontSize',24,'Location','northeast','AutoUpdate','off');
% psfrTools.makeAxisSquared(h);
% plot([xlim()],[1,1],'k:','LineWidth',1.2)
% text(5,1.1,'1\%-FVU','interpreter','latex','FontSize',24);

%% ---------------------- Section 4.2: FOV = f(h) ---------------------- %%

nL= 100;
close all;
% Atmosphere defintion
if nL == 7
    zenithAngle = 45;
    airmass     = 1/cos(zenithAngle*pi/180);
    r0          = 0.157*airmass^(-3/5);
    h0          = [100, 500, 1000, 2000, 4000, 8000, 16000]*airmass;
    w0          = r0^(-5/3)*[60.4, 11.3, 9.6, 2.9, 5.8, 4.2, 5.8]/100;
    L0          = 25;
    wS          = ones(1,nL);
    wD          = zeros(1,nL);  
elseif nL == 35
    
    profile = fitsread('profil_turbulent_eso.fits');
    h0     = profile(1:35,2)';
    w0     = profile(1:35,4)/100;
    w0     = r0^(-5/3)*w0/sum(w0);
    wS    = profile(1:35,3);
    wD    = profile(1:35,9);
else
    zenithAngle = 45;
    airmass     = 1/cos(zenithAngle*pi/180);
    r0          = 0.157*airmass^(-3/5);
    L0          = 25;
    wS          = ones(1,nL);
    wD          = zeros(1,nL); 
    w0          = r0^(-5/3)*ones(1,nL)/nL;
    h0          = logspace(log10(30),log10(20e3),nL);  
end

atm = atmosphere(photometry.V0,r0,L0,'altitude',h0,'fractionnalR0',w0/sum(w0),...
    'windSpeed',wS,'windDirection',wD);
    
% Inputs parameters
tel.obstructionRatio = 0.3;
Dvec = [0.5,1,2,4];
nD   = length(Dvec);
%wvl  = (400:50:900)*1e-9;
%nWvl = length(wvl);
redo = false;

if redo        
    aspRatio = zeros(nD,nL,nSrc);
    FWHM     = zeros(nD,nL,nSrc);
    FOV      = zeros(nD,nL);
    nSrc     = 100;
    aS       = pi/4*ones(1,nSrc);
    pp = [1.5821    2.5104    6.6322   10.2994];

    tic;
    for d=1:nD           % Telescope diameter 
        tel.D = Dvec(d);   
        for l=1:nL   % Bin height
                % Update the FOV to reduce the number of sources
                zS    = constants.arcmin2radian*pp(d)*(h0(l)/1e3)^(-1)*linspace(0.995,1.005,nSrc);                
                src2  = source('wavelength',photometry.V0,'zenith',zS,'azimuth',aS);                
                ps    = constants.radian2mas*src2(1).wavelength/tel.D/2;
                %Instantiate the FPP
                fppr = focalPlaneProfiler(tel,src2,ngs(1),atm.slab(l),otf0,ps,otf0,...
                    [],[],'nLayer',1,'weightInit',w0(l),'fitWeights',true,'heightInit',h0(l),...
                    'fitHeight',false,'nptZonal',dk,'idxPhase',validPhase,'flagAnisokinetism',true);
                
                % Calculate the PSF
                psfMod = fppr.psfModel(w0(l),fppr.otfKernel);
                % Measure the aspect ratio
                for iSrc = 1:nSrc
                    idSrc = 1 + (iSrc-1)*npsf:iSrc*npsf;
                    [FWHMx,FWHMy,aRatio] = psfrTools.getFWHM(psfMod(:,idSrc),ps,2);
                    FWHM(d,l,iSrc)     = max([FWHMx,FWHMy]);
                    aspRatio(d,l,iSrc) = aRatio;                    
                end
                % Estimate the FOV regarding the maximal PSF elongation
                FOV(d,l) = constants.radian2arcmin*zS(find(aspRatio(d,l,:) == max(aspRatio(d,l,:))));
        end
    end
    t0 = toc();
    
    fitswrite(aspRatio,[savePath,'fits/psfElipticityvDistvAltvDtelvWavelength_',int2str(nL),'layers.fits']);
    fitswrite(FWHM,[savePath,'fits/psfFWHMvDistvAltvDtelvWavelength_',int2str(nL),'layers.fits']);
    fitswrite(FOV,[savePath,'fits/FOVvAltvDtelvWavelength_',int2str(nL),'layers.fits']);
else
    aspRatio = fitsread([savePath,'fits/psfElipticityvDistvAltvDtelvWavelength_',int2str(nL),'layers.fits']);
    FWHM     = fitsread([savePath,'fits/psfFWHMvDistvAltvDtelvWavelength_',int2str(nL),'layers.fits']);
    FOV      = fitsread([savePath,'fits/FOVvAltvDtelvWavelength_',int2str(nL),'layers.fits']);
end


% Plot theta = f(h)
close all;
x = h0/1e3;
F = @(a,x) a*x;
p = zeros(1,nD);
dp=zeros(1,nD);
for d=1:nD
     [p(d),~,R,~,~,~,J] = lsqcurvefit(F,1,x.^(-1),FOV(d,:));
     dp(d) = diff(nlparci(p(d),R,'jacobian',J),1,2);
end
pa = mean(p./Dvec);
dpa= sqrt(mean((dp(d)./Dvec).^2)) + std(p./Dvec)/sqrt(nD);
idx = 1:5:length(x);

h = figure;
loglog(x(idx),FOV(1,idx),'ks','MarkerFaceColor','k','MarkerSize',7,'LineWidth',1.2);hold on;
loglog(x(idx),FOV(2,idx),'bo','MarkerFaceColor','b','MarkerSize',7,'LineWidth',1.2);
loglog(x(idx),FOV(3,idx),'rd','MarkerFaceColor','r','MarkerSize',7,'LineWidth',1.2);
loglog(x(idx),FOV(4,idx),'mp','MarkerFaceColor','m','MarkerSize',7,'LineWidth',1.2);
xlabel('Altitude (km)','interpreter','latex','FontSize',24);
ylabel('FOV (arcmin)','interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
legend({'50cm-telescope','1m-telescope','2m-telescope','4m-telescope'},...
    'interpreter','latex','FontSize',24,'Location','northeast','AutoUpdate','off');
% Plot the linear regression
plot(x,p(1)*x.^(-1),'k--');
plot(x,p(2)*x.^(-1),'b--');
plot(x,p(3)*x.^(-1),'r--');
plot(x,p(4)*x.^(-1),'m--');
psfrTools.makeAxisSquared(h)
xlim([0.95*min(x),max(x)*1.05])
ylim([0.95*min(FOV(:)),max(FOV(:))*1.05])
text(1.2*min(x),1.2*min(FOV(:)),['FOV[arcmin] =$',num2str(pa),'\times D[m]/h[km]$']...
    ,'interpreter','latex','FontSize',24);


%% ---------------------- Section 4.3: dh = f(h) ---------------------- %%
close all;
% Altitude resolution versus height
dt1 = (0.1/60);
dt2 = (0.5/60);
dt3 = (1/60);
hh  = linspace(0.1,20,1000);
dh1 = dt1/p(1)*hh.^2;
dh2 = dt2/p(1)*hh.^2;
dh3 = dt3/p(1)*hh.^2;

h = figure;
loglog(hh,dh1*1e3,'k--','LineWidth',1.2);hold on
loglog(hh,dh2*1e3,'b-.','LineWidth',1.2)
loglog(hh,dh3*1e3,'r-','LineWidth',1.2)
legend({'$\Delta\theta=0.1"$','$\Delta\theta=0.5"$','$\Delta\theta=1"$'},'interpreter','latex','FontSize',24,...
    'Location','north','AutoUpdate','off');
psfrTools.makeAxisSquared(h)
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
xlabel('Altitude (km)','interpreter','latex','FontSize',24);
ylabel('Altitude resolution (m) for D=0.5m','interpreter','latex','FontSize',24);
ylim([min(dh1*1e3),max(dh3*1e3)]);
xlim([min(hh),max(hh)]);
plot([0.2,0.2],[ylim()],'k:');
plot([4,4],[ylim()],'k:');
plot([12,12],[ylim()],'k:');
text(0.15,1.5e-2,'100m','interpreter','latex','FontSize',24);
text(3,1.5e-2,'4km','interpreter','latex','FontSize',24);
text(10,1.5e-2,'12km','interpreter','latex','FontSize',24);
plot([xlim()],[1,1],'k:');
plot([xlim()],[100,100],'k:');
plot([xlim()],[1e3,1e3],'k:');
text(1.1e-1,1.3,'1m','interpreter','latex','FontSize',24);
text(1.1e-1,130,'100m','interpreter','latex','FontSize',24);
text(1.1e-1,1300,'1km','interpreter','latex','FontSize',24);

%% -------------------- Section 4.4: Layers completness ---------------- %%
redo =false;
close all;
if redo
    %atmosphere
    nL   = 7;
    zenithAngle = 45;
    airmass     = 1/cos(zenithAngle*pi/180);
    r0          = 0.157*airmass^(-3/5);
    h0          = [100, 500, 1000, 2000, 4000, 8000, 16000]*airmass;
    w0          = r0^(-5/3)*[60.4, 11.3, 9.6, 2.9, 5.8, 4.2, 5.8]/100;
    L0          = 25;
    wS          = ones(1,nL);
    wD          = zeros(1,nL);
    % Sources
    nSrc = 480*4+1;
    zS   = linspace(0,480,nSrc)*constants.arcsec2radian;
    aS   = pi/4*ones(1,nSrc);
    src  = source('wavelength',photometry.V0,'zenith',zS,'azimuth',aS);
    % Get the optimal position regarding the maximal PSF elongation
    aspRatio = fitsread([savePath,'fits/psfElipticityvDistvAlt_7layers.fits']);
    idM = zeros(1,nL);
    for l=1:nL
        idM(l) = find(aspRatio(l,:) == max(aspRatio(l,:)));
    end
    
    nK  = 100;
    Wk  = zeros(nL,nK);
    fW  = 0*Wk;
    
    for k=1:nK
        fppr = focalPlaneProfiler(tel,src(idM+k-1),ngs(1),atm,otf0,psInMas,otf0,...
            [],[],'nLayer',nL,'weightInit',w0,'fitWeights',true,'heightInit',h0,...
            'fitHeight',false,'nptZonal',dk,'idxPhase',validPhase,'flagAnisokinetism',true);
        
        psfModel  = fppr.psfModel(w0,fppr.otfKernel);
        
        % Cn2h retrieval
        fppr.psfSky     = psfModel;
        fppr.weightInit = ones(1,nL)/nL*atm.r0^(-5/3);
        fppr            = fppr.getParameters();
        Wk(:,k)        = fppr.xfinal;
        fW(:,k)        = fppr.dweights;
    end
    dWk = 100*bsxfun(@rdivide,abs(bsxfun(@minus,Wk,w0')),w0');
    fitswrite([dWk,fW],[savePath,'fits/cn2hErrorVshift_',int2str(nL),'layers.fits']);
else
    out = fitsread([savePath,'fits/cn2hErrorVshift_',int2str(nL),'layers.fits']);
    dWk = out(:,1:end/2);
    fW  = out(:,end/2+1:end);
    [nL,nK]= size(dWk);
end

% Plot cn2h error =f(delta theta)
th_anisoTT  = 0.5*(zernikeStats.anisokinetismAngle(zernike(2),atm) + zernikeStats.anisokinetismAngle(zernike(3),atm));

% Define altitude range
idLow  = h0<1e3;
idMid  = (h0>=1e3) & (h0<8e3);
idHigh = (h0>=8e3) & (h0<16e3);
idVhigh= h0>=16e3;


dT= mean(diff(zS))*constants.radian2arcsec*(0:nK-1);
nn = 50;
%Low altitude bins <1km
dW1= smooth(sqrt(mean(dWk(idLow,:).^2,1)),nn)';
eW1= smooth(3*sqrt(mean(fW(idLow,:).^2,1)),nn)';
%Mid altitude bins <8km
dW2= smooth(sqrt(mean(dWk(idMid,:).^2,1)),nn)';
eW2= smooth(3*sqrt(mean(fW(idMid,:).^2,1)),nn)';
%High altitude bins <16km
dW3= smooth(sqrt(mean(dWk(idHigh,:).^2,1)),nn)';
eW3= smooth(3*sqrt(mean(fW(idHigh,:).^2,1)),nn)';
%High altitude bins <16km
dW4= smooth(sqrt(mean(dWk(idVhigh,:).^2,1)),nn)';
eW4= smooth(3*sqrt(mean(fW(idVhigh,:).^2,1)),nn)';


h = figure;
semilogy(dT,dW1,'ks--','MarkerFaceColor','k','MarkerSize',3,'LineWidth',1.2);hold on;
semilogy(dT,dW2,'bo--','MarkerFaceColor','b','MarkerSize',3,'LineWidth',1.2);
semilogy(dT,dW3,'rd--','MarkerFaceColor','r','MarkerSize',3,'LineWidth',1.2);
semilogy(dT,dW4,'mp--','MarkerFaceColor','m','MarkerSize',3,'LineWidth',1.2);
xlabel('Positioning shift (arcsec)','interpreter','latex','FontSize',24);
ylabel('$C_n^2(h)$ error (\%)','interpreter','latex','FontSize',24);
legend({'$h< 1$km','8km$> h\geq$1km','16km$>h\geq$8km','$h\geq 16$km'},'interpreter','latex',...
    'FontSize',24,'Location','southeast','AutoUpdate','off');
psfrTools.makeAxisSquared(h);
inBetween = [(dW1-eW1) fliplr(dW1+eW1)];
fill([dT fliplr(dT)]', inBetween, 'k','FaceAlpha',.15,'LineStyle','none');
inBetween = [(dW2-eW2) fliplr(dW2+eW2)];
fill([dT fliplr(dT)]', inBetween, 'b','FaceAlpha',.15,'LineStyle','none');
inBetween = [(dW3-eW3) fliplr(dW3+eW3)];
fill([dT fliplr(dT)]', inBetween, 'r','FaceAlpha',.15,'LineStyle','none');
inBetween = [(dW4-eW4) fliplr(dW4+eW4)];
fill([dT fliplr(dT)]', inBetween, 'm','FaceAlpha',.15,'LineStyle','none');
xlim([0,max(dT)])
ylim([0.001,20])
loglog([xlim()],[1,1],'k--');
text(1,1.2,'1\%','interpreter','latex','FontSize',24);
plot([1,1]*th_anisoTT,[ylim()],'k--');
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
text(1.02*th_anisoTT,2e-3,'$\theta_a$','interpreter','latex','FontSize',24);


%% ------------------- Section 4.5: Noise = f(D,psInMas) - ------------- %%
close all;
redo =false;

fov = 1;
nRes= 100;
Fe  = 1000;
tel = telescope(0.5,'fieldOfViewInArcMin',fov,'resolution',nRes,...
    'samplingTime',1/Fe,'obstructionRatio',0.3);
ngs0 =  source('wavelength',photometry.V0);
h1   = 10e3;
w1   = 30.9515;
atm2 = atmosphere(photometry.V0,w1^(-3/5),25,'altitude',h1,'fractionnalR0',1,...
    'windSpeed',0,'windDirection',0);
src2 = source('wavelength',photometry.V0,'zenith',12*constants.arcsec2radian,'azimuth',pi/4);
tel.obstructionRatio = 0.3;
% CCD config
tExp      = 1e-2;
nExp      = 3000;
ron       = 0.8;
darkBg    = 0.05;
eCCD      = sqrt(nExp*(ron^2 + (darkBg*tExp)^2));
QE        = 0.7;
throughput= 0.9;
deltaL    = 0.1*300;%300nm -> amgstrom
F0        = 995.5*1e6*deltaL;%3.3e12;
%995.5 ph/cm^2/s
mag       = 10:16;
nMag      = length(mag);
nIter     = 100;
Dvec      = [0.5,1,2,4];
nD        = length(Dvec);
nq        = [2,1,0.5,0.4,0.3,0.25,0.24,0.23,0.22,0.21,0.2,0.19,0.18,0.15,0.1];
nP        = length(nq);
Wn        = zeros(nD,nP,nMag,nIter);
fW        = 0*Wn;
        
if redo 
    tic;
    for d=1:nD
        % Changing the telescope diameter
        tel.D = Dvec(d);
        F2DN  = tExp*nExp*throughput*QE*tel.area;
        fov   = 25*constants.radian2mas*src2.wavelength/tel.D;
        for p=1:nP
            % Changing the camera pixel scale
            ps       = constants.radian2mas*src2.wavelength/tel.D/nq(p)/2;
            skyBg    = (ps*1e-3)^2*10^(-0.4*24)*F0*F2DN;
            fovInPix = round(fov/ps)+1;
            % Instantiate the FPPR
            fppr   = focalPlaneProfiler(tel,src2,ngs0,atm2,zeros(fovInPix),ps,otf0,...
                [],[],'nLayer',1,'weightInit',w1,'fitWeights',true,'heightInit',h1,...
                'fitHeight',false,'nptZonal',dk,'idxPhase',validPhase,...
                'flagAnisokinetism',true);
            
            psfModel  = fppr.psfModel(w1,fppr.otfKernel);
            psfModel  = psfModel/sum(psfModel(:));
                                              
            for k=1:nMag
                for t=1:nIter
                    % Noisify the image                    
                    tmp = psfModel*10^(-0.4*mag(k))*F0*F2DN;
                    fppr.psfSky = poissrnd(tmp + skyBg) - skyBg + eCCD*randn(fovInPix);
                    
                    % Cn2h retrieval
                    fppr.weightInit = 15;
                    fppr            = fppr.getParameters();
                    Wn(d,p,k,t)     = fppr.xfinal;
                    fW(d,p,k,t)     = fppr.dweights;
                end
            end
        end
    end
    t0=toc;
    dWn = 100*(Wn-w1)/w1;    
    fitswrite([dWn,fW],[savePath,'fits/cn2hError_Dtel_pscale_Mag_1layer_fovconstant_largeband.fits']);        
else
    out = fitsread([savePath,'fits/cn2hError_Dtel_pscale_Mag_1layer_fovconstant.fits']);
    dWn = out(:,1:end/2,:,:);
    fW  = out(:,end/2+1:end,:,:);
    [nD,nP,nMag,nIter] = size(dWn);
    mag = 8:1:14;    
    %mag = 10:16;
end



%% Plot Cn2h error v magnitude
close all;

dW = mean(abs(dWn),4);
eW = sqrt(var(dWn,[],4) + mean(fW.^2,4))/sqrt(nIter);
% error = f(mag,D);
dW1= squeeze(dW(1,11,:))';
eW1= squeeze(eW(1,11,:))';
dW2= squeeze(dW(2,11,:))';
eW2= squeeze(eW(2,11,:))';
dW3= squeeze(dW(3,11,:))';
eW3= squeeze(eW(3,11,:))';
dW4= squeeze(dW(4,11,:))';
eW4= squeeze(eW(4,11,:))';

h = figure;
plot(mag,dW1,'ks--','MarkerFaceColor','k','MarkerSize',7,'LineWidth',1.2);hold on;
plot(mag,dW2,'bo--','MarkerFaceColor','b','MarkerSize',7,'LineWidth',1.2);
plot(mag,dW3,'rd--','MarkerFaceColor','r','MarkerSize',7,'LineWidth',1.2);
plot(mag,dW4,'mp--','MarkerFaceColor','m','MarkerSize',7,'LineWidth',1.2);
xlabel('Calibration star magnitude (mag)','interpreter','latex','FontSize',24);
ylabel('$C_n^2(h)$ error (\%)','interpreter','latex','FontSize',24);
legend({'50cm-telescope','1m-telescope','2m-telescope','4m-telescope'},...
    'interpreter','latex','FontSize',24,'Location','northwest','AutoUpdate','off');
psfrTools.makeAxisSquared(h);
inBetween = [(dW1-eW1) fliplr(dW1+eW1)];
fill([mag fliplr(mag)], inBetween, 'k','FaceAlpha',.15,'LineStyle','none');
inBetween = [(dW2-eW2) fliplr(dW2+eW2)];
fill([mag fliplr(mag)], inBetween, 'b','FaceAlpha',.15,'LineStyle','none');
inBetween = [(dW3-eW3) fliplr(dW3+eW3)];
fill([mag fliplr(mag)], inBetween, 'r','FaceAlpha',.15,'LineStyle','none');
inBetween = [(dW4-eW4) fliplr(dW4+eW4)];
fill([mag fliplr(mag)], inBetween, 'm','FaceAlpha',.15,'LineStyle','none');
xlim([mag(1)-1,mag(end)+1]);
plot([xlim()],[1,1],'k--');
text(7.5,1.1,'1\%','interpreter','latex','FontSize',24);
%plot([xlim()],[10,10],'k--');
%text(7.5,11.8,'10\%','interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
set(gca,'YScale','log')

% error = f(mag,ps);
dW1= squeeze(dW(1,:,4));
eW1= squeeze(eW(1,:,4));
dW2= squeeze(dW(1,:,5));
eW2= squeeze(eW(1,:,5));
dW3= squeeze(dW(1,:,6));
eW3= squeeze(eW(1,:,6));
dW4= squeeze(dW(1,:,7));
eW4= squeeze(eW(1,:,7));

pp = constants.radian2mas*src2.wavelength/0.5./nq/2;
h = figure;
plot(pp,dW1,'ks--','MarkerFaceColor','k','MarkerSize',7,'LineWidth',1.2);hold on;
plot(pp,dW2,'bo--','MarkerFaceColor','b','MarkerSize',7,'LineWidth',1.2);
plot(pp,dW3,'rd--','MarkerFaceColor','r','MarkerSize',7,'LineWidth',1.2);
plot(pp,dW4,'mp--','MarkerFaceColor','m','MarkerSize',7,'LineWidth',1.2);
xlabel('Pixel scale (mas)','interpreter','latex','FontSize',24);
ylabel('$C_n^2(h)$ error (\%)','interpreter','latex','FontSize',24);
legend({'V=13','V=14','V=15','V=16'},...
    'interpreter','latex','FontSize',24,'Location','northeast','AutoUpdate','off');
psfrTools.makeAxisSquared(h);
inBetween = [(dW1-eW1) fliplr(dW1+eW1)];
fill([pp fliplr(pp)], inBetween, 'k','FaceAlpha',.15,'LineStyle','none');
inBetween = [(dW2-eW2) fliplr(dW2+eW2)];
fill([pp fliplr(pp)], inBetween, 'b','FaceAlpha',.15,'LineStyle','none');
inBetween = [(dW3-eW3) fliplr(dW3+eW3)];
fill([pp fliplr(pp)], inBetween, 'r','FaceAlpha',.15,'LineStyle','none');
inBetween = [(dW4-eW4) fliplr(dW4+eW4)];
fill([pp fliplr(pp)], inBetween, 'm','FaceAlpha',.15,'LineStyle','none');
xlim([20,1050]);
plot([xlim()],[1,1],'k--');
text(30,1.05,'1\%','interpreter','latex','FontSize',24);
plot([xlim()],[5,5],'k--');
text(30,5.3,'5\%','interpreter','latex','FontSize',24);
set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
set(gca,'YScale','log')


%% OLD

%% FOV = f(theta,D,lambda) model
% Dvec = [0.5,1,2,4];
% nD   = length(Dvec);
% wvl  = (400:50:900)*1e-9;
% nWvl = length(wvl);
% nSrc2= 101;
% zS   = linspace(0,100,nSrc2)*constants.arcsec2radian;
% aS   = pi/4*ones(1,nSrc2);
% h1   = 8e3;
% atm2 = atmosphere(photometry.V0,r0,L0,'altitude',h1,'fractionnalR0',1,...
%     'windSpeed',0,'windDirection',0);
% 
% aspRatio = zeros(nD,nWvl,nSrc2);
% FWHM     = zeros(nD,nWvl,nSrc2);
% maxE     = zeros(nD,nWvl);
% 
% redo  =false;
% 
% if redo
%     for d=1:nD
%         tel.D = Dvec(d);   
%         tel.obstructionRatio = 0.3;
%         for l=1:nWvl
%             % Define a new source
%             src2  = source('wavelength',photometryKASP(wvl(l),0,0),'zenith',zS,'azimuth',aS);
%             atm2.wavelength = src2(1).wavelength;
%             w1    = atm2.r0^(-5/3);
%             % Define the pixel scale
%             ps    = constants.radian2mas*wvl(l)/tel.D/4;
%             %Instantiate the FPPR
%             fppr  = focalPlaneProfiler(tel,src2,ngs(1),atm2,otf0,ps,otf0,...
%                 [],[],'nLayer',atm2.nLayer,'weightInit',w1,'fitWeights',true,'heightInit',h1,...
%                 'fitHeight',false,'nptZonal',dk,'idxPhase',validPhase,'flagAnisokinetism',true);
%             % Create the PSF
%             psfMod = fppr.psfModel(w1,fppr.otfKernel);
%             % Evaluate the PSF aspect ratio/FWHM
%             for iSrc = 1:nSrc2
%                 idSrc = 1 + (iSrc-1)*npsf:iSrc*npsf;
%                 [FWHMx,FWHMy,aRatio] = psfrTools.getFWHM(psfMod(:,idSrc),ps,2);
%                 FWHM(d,l,iSrc)         = max([FWHMx,FWHMy]);
%                 aspRatio(d,l,iSrc)     = aRatio;
%             end
%             maxE(d,l) = zS(find(aspRatio(d,l,:) == max(aspRatio(d,l,:))))*constants.radian2arcmin;
%         end        
%     end
%     fitswrite(aspRatio,[savePath,'fits/aspectRatiovDvLambda_1layer.fits']);
%     fitswrite(FWHM,[savePath,'fits/fwhmDvLambda_1layer.fits']);
%     fitswrite(maxE,[savePath,'fits/argmaxAspectRatiovDvLambda_1layer.fits']);
% 
% else
%     aspRatio = fitsread([savePath,'fits/aspectRatiovDvLambda_1layer.fits']);
%     FWHM     = fitsread([savePath,'fits/fwhmDvLambda_1layer.fits']);
%     maxE     = fitsread([savePath,'fits/argmaxAspectRatiovDvLambda_1layer.fits']);
% end
% 
% % Aspect ratio = f(D,theta)
% h = figure;
% semilogx(zS*constants.radian2arcsec,squeeze(aspRatio(1,3,:)),'k--','LineWidth',1.2);
% hold on;
% semilogx(zS*constants.radian2arcsec,squeeze(aspRatio(2,3,:)),'b:','LineWidth',1.2);
% semilogx(zS*constants.radian2arcsec,squeeze(aspRatio(3,3,:)),'r-.','LineWidth',1.2);
% semilogx(zS*constants.radian2arcsec,squeeze(aspRatio(4,3,:)),'m-','LineWidth',1.2);
% ylabel('PSF Aspect ratio in V-band','interpreter','latex','FontSize',24);
% xlabel('Off-axis distance (arcsec)','interpreter','latex','FontSize',24);
% set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
% legend({'50cm-telescope','1m-telescope','2m-telescope','4m-telescope'},...
%     'interpreter','latex','FontSize',24,'Location','northwest');
% psfrTools.makeAxisSquared(h);
% 
% % Aspect ratio = f(wvl,theta) 4-m telescope
% h = figure;
% semilogx(zS*constants.radian2arcsec,squeeze(aspRatio(4,2,:)),'k--','LineWidth',1.2);
% hold on;
% semilogx(zS*constants.radian2arcsec,squeeze(aspRatio(4,3,:)),'b:','LineWidth',1.2);
% semilogx(zS*constants.radian2arcsec,squeeze(aspRatio(4,6,:)),'r-.','LineWidth',1.2);
% semilogx(zS*constants.radian2arcsec,squeeze(aspRatio(4,10,:)),'m-','LineWidth',1.2);
% ylabel('PSF Aspect ratio on a 4m-class telescope','interpreter','latex','FontSize',24);
% xlabel('Off-axis distance (arcsec)','interpreter','latex','FontSize',24);
% set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
% legend({'B-band','V-band','R-band','I-band'},...
%     'interpreter','latex','FontSize',24,'Location','northwest');
% psfrTools.makeAxisSquared(h);
% 
% 
% % Aspect ratio = f(wvl,theta)
% h = figure;
% semilogx(wvl*1e9,squeeze(maxE(1,:))*60,'ks--','MarkerFaceColor','k','MarkerSize',5,'LineWidth',1.2);
% hold on;
% semilogx(wvl*1e9,squeeze(maxE(2,:))*60,'bo:','MarkerFaceColor','b','MarkerSize',5,'LineWidth',1.2);
% semilogx(wvl*1e9,squeeze(maxE(3,:))*60,'rd-.','MarkerFaceColor','r','MarkerSize',5,'LineWidth',1.2);
% semilogx(wvl*1e9,squeeze(maxE(4,:))*60,'mp-','MarkerFaceColor','m','MarkerSize',5,'LineWidth',1.2);
% xlabel('Wavelength (nm)','interpreter','latex','FontSize',24);
% ylabel('FOV (arcmin)','interpreter','latex','FontSize',24);
% set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
% legend({'50cm-telescope','1m-telescope','2m-telescope','4m-telescope'},...
%     'interpreter','latex','FontSize',24,'Location','northwest');
% psfrTools.makeAxisSquared(h);



% if nL == 7
%     h = figure;
%     semilogx(zS*constants.radian2arcmin,aspRatio(3,:),'k:','LineWidth',1.2);hold on;   
%     semilogx(zS*constants.radian2arcmin,aspRatio(5,:),'b--','LineWidth',1.2);
%     semilogx(zS*constants.radian2arcmin,aspRatio(6,:),'r-.','LineWidth',1.2);
%     semilogx(zS*constants.radian2arcmin,aspRatio(7,:),'m-','LineWidth',1.2);
%     ylabel('Aspect ratio','interpreter','latex','FontSize',24);
%     xlabel('Off-axis distance (arcmin)','interpreter','latex','FontSize',24);
%     set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
%     legend({'h=1km','h=4km','h=8km','h=16km'},...
%         'interpreter','latex','FontSize',24,'Location','northwest');
%     psfrTools.makeAxisSquared(h)
%     [m,n] = psfrTools.minmax(zS*constants.radian2arcmin);
%     xlim([m,n])
%     
%     h = figure;
%     semilogx(zS*constants.radian2arcmin,FWHM(3,:),'k:','LineWidth',1.2);hold on;   
%     semilogx(zS*constants.radian2arcmin,FWHM(5,:),'b--','LineWidth',1.2);
%     semilogx(zS*constants.radian2arcmin,FWHM(6,:),'r-.','LineWidth',1.2);
%     semilogx(zS*constants.radian2arcmin,FWHM(7,:),'m-','LineWidth',1.2);
%     ylabel('FWHM (mas)','interpreter','latex','FontSize',24);
%     xlabel('Off-axis distance (arcmin)','interpreter','latex','FontSize',24);
%     set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
%     legend({'h=1km','h=4km','h=8km','h=16km'},...
%         'interpreter','latex','FontSize',24,'Location','northwest');
%     psfrTools.makeAxisSquared(h)
%     xlim([m,n])
%     
%     idx = [2:7];
% else
%     h = figure;
%     semilogx(zS*constants.radian2arcmin,aspRatio(8,:));hold on;
%     semilogx(zS*constants.radian2arcmin,aspRatio(9,:));
%     semilogx(zS*constants.radian2arcmin,aspRatio(11,:));
%     semilogx(zS*constants.radian2arcmin,aspRatio(14,:));
%     semilogx(zS*constants.radian2arcmin,aspRatio(19,:));
%     semilogx(zS*constants.radian2arcmin,aspRatio(25,:));
%     semilogx(zS*constants.radian2arcmin,aspRatio(29,:));
%     ylabel('Aspect ratio','interpreter','latex','FontSize',24);
%     xlabel('Off-axis distance (arcmin)','interpreter','latex','FontSize',24);
%     set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
%     legend({'h=0.6km','h=1.13km','h=2.63km','h=5.5km','h=10.5km','h=16.5km','h=20.5km'},...
%         'interpreter','latex','FontSize',24,'Location','southwest');
%     psfrTools.makeAxisSquared(h)
%     
%     
%     
%     h = figure;
%     semilogx(zS*constants.radian2arcmin,FWHM(8,:));hold on;
%     semilogx(zS*constants.radian2arcmin,FWHM(9,:));
%     semilogx(zS*constants.radian2arcmin,FWHM(11,:));
%     semilogx(zS*constants.radian2arcmin,FWHM(14,:));
%     semilogx(zS*constants.radian2arcmin,FWHM(19,:));
%     semilogx(zS*constants.radian2arcmin,FWHM(25,:));
%     semilogx(zS*constants.radian2arcmin,FWHM(29,:));
%     ylabel('FWHM (mas)','interpreter','latex','FontSize',24);
%     xlabel('Off-axis distance (arcmin)','interpreter','latex','FontSize',24);
%     set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex' );
%     legend({'h=0.6km','h=1.13km','h=2.63km','h=5.5km','h=10.5km','h=16.5km','h=20.5km'},...
%         'interpreter','latex','FontSize',24,'Location','northwest');
%     psfrTools.makeAxisSquared(h)
%     
%     idx = [5:35];
% end


% 
% close all
% redo = false;
% if redo    
%     Wl       = zeros(nL,nL);
%     fW       = 0*Wl;
%     
%     for l=1:nL        
%         msk = idM(l:end);
%         fppr = focalPlaneProfiler(tel,src(msk),ngs(1),atm,otf0,psInMas,otf0,...
%             [],[],'nLayer',nL,'weightInit',w0,'fitWeights',true,'heightInit',h0,...
%             'fitHeight',false,'nptZonal',dk,'idxPhase',validPhase,'flagAnisokinetism',true);
%         psfModel  = fppr.psfModel(w0,fppr.otfKernel);
%         
%         % Cn2h retrieval
%         fppr.psfSky     = psfModel;
%         fppr.weightInit = ones(1,nL)/nL*atm.r0^(-5/3);
%         fppr            = fppr.getParameters();
%         Wl(:,l)        = fppr.xfinal;
%         fW(:,l)        = fppr.dweights;
%     end
%     dWl = 100*bsxfun(@rdivide,abs(bsxfun(@minus,Wl,w0')),w0');
%     fitswrite([dWl,fW],[savePath,'fits/cn2hErrorVfov_',int2str(nL),'layers.fits']);
% else
%     out = fitsread([savePath,'fits/cn2hErrorVfov_',int2str(nL),'layers.fits']);
%     dWl = out(:,1:end/2);
%     fW  = out(:,end/2+1:end);
%     nL  = size(dWl,1);
%     
% end
% 
% % Plot Cn2h error =f(FOV)
% %fovArray = linspace(480,3,nFov);
% 
% fovArray= zS(idM)*constants.radian2arcsec/60;
% nn = 50;
% %Low altitude bins <1km
% dW1= smooth(sqrt(mean(dWl(idLow,:).^2,1)),nn)';
% %Mid altitude bins <8km
% dW2= smooth(sqrt(mean(dWl(idMid,:).^2,1)),nn)';
% %High altitude bins <16km
% dW3= smooth(sqrt(mean(dWl(idHigh,:).^2,1)),nn)';
% %High altitude bins <16km
% dW4= smooth(sqrt(mean(dWl(idVhigh,:).^2,1)),nn)';
% 
% 
% h = figure;
% loglog(fovArray,dW1,'ks--','MarkerFaceColor','k','MarkerSize',7,'LineWidth',1.2);hold on;
% loglog(fovArray,dW2,'bo--','MarkerFaceColor','b','MarkerSize',7,'LineWidth',1.2);
% loglog(fovArray,dW3,'rd--','MarkerFaceColor','r','MarkerSize',7,'LineWidth',1.2);
% loglog(fovArray,dW4,'mp--','MarkerFaceColor','m','MarkerSize',7,'LineWidth',1.2);
% xlabel('Field of view (arcmin)','interpreter','latex','FontSize',24);
% ylabel('$C_n^2(h)$ error (\%)','interpreter','latex','FontSize',24);
% legend({'$h< 1$km','8km$> h\geq$1km','16km$>h\geq$8km','$h\geq 16$km'},'interpreter','latex',...
%     'FontSize',24,'Location','southwest','AutoUpdate','off');
% psfrTools.makeAxisSquared(h);
% plot([xlim()],[1,1],'k--');
% text(1,1.2,'1\%','interpreter','latex','FontSize',24);
% plot([xlim()],[10,10],'k--');
% text(1,11.5,'10\%','interpreter','latex','FontSize',24);
% set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
% 


%%Section 3.2: Cn2h Retrieval performace
% redo = false;
% 
% if redo
%     nT = 100;
%     Wt = zeros(nL,nT);
%     fW = 0*Wt;
%     
%     fppr = focalPlaneProfiler(tel,src(idM),ngs(1),atm,otf0,psInMas,otf0,...
%         [],[],'nLayer',nL,'weightInit',w0,'fitWeights',true,'heightInit',h0,...
%         'fitHeight',false,'nptZonal',dk,'idxPhase',validPhase,'flagAnisokinetism',true);
%     psfModel       = fppr.psfModel(w0,fppr.otfKernel);
%     fppr.psfSky     = psfModel;
%             
%     for t=1:nT
%         fppr.weightInit = atm.r0^(-5/3).*abs(randn(1,nL));
%         fppr            = fppr.getParameters();
%         Wt(:,t)        = fppr.xfinal;
%         fW(:,t)        = fppr.dweights;
%     end
%     fitswrite([Wt,fW],[savePath,'fits/cn2hErrorVinitGuess_',int2str(nL),'layers.fits']);
% else
%     out = fitsread([savePath,'fits/cn2hErrorVinitGuess_',int2str(nL),'layers.fits']);
%     Wt  = out(:,1:end/2);
%     fW  = out(:,end/2+1:end);
%     [nL,nT]= size(Wt);
% end
% 
% 
% 
% dWt = 1e2*bsxfun(@rdivide,(bsxfun(@minus,Wt,w0')),w0');
% h = figure;
% histogram(dWt(:),30,'Normalization','probability');
% ylabel('Probability','interpreter','latex','FontSize',24);
% xlabel('$C_n^2$ error (\%))','interpreter','latex','FontSize',24);
% psfrTools.makeAxisSquared(h);
% set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');

%mean:0.017\%
%std: 0.18\% 1-sigma

% %% ------------------- Section 4.3: Noise = f(Mag) --------------------- %%
% redo =false;
% 
% if redo
%     fppr = focalPlaneProfiler(tel,src(idM),ngs(1),atm,otf0,psInMas,otf0,...
%         [],[],'nLayer',nL,'weightInit',w0,'fitWeights',true,'heightInit',h0,...
%         'fitHeight',false,'nptZonal',dk,'idxPhase',validPhase,'flagAnisokinetism',true);
%     
%     psfModel  = fppr.psfModel(w0,fppr.otfKernel);
%     psfModel  = psfModel/sum(psfModel(:));
%     psfNoise  = 0*psfModel;
%     
%     % CCD config
%     tExp      = 1e-2;
%     nExp      = 3000;
%     ron       = 0.8;
%     darkBg    = 0.05;
%     eCCD      = sqrt(nExp*(ron^2 + (darkBg*tExp)^2));
%     QE        = 0.7;
%     throughput= 0.9;
%     F0        = 3.3e12;
%     %995.5 ph/cm^2/s
%     F2DN      = tExp*nExp*throughput*QE*tel.area;
%     skyBg     = (psInMas*1e-3)^2*10^(-0.4*24)*F0*F2DN;
%     mag       = 0:13;
%     nMag      = length(mag);
%     nIter     = 100;
%     Wn        = zeros(nL,nMag,nIter);
%     fW        = 0*Wn;
%     
%     tic;
%     for k=1:nMag
%         for t=1:nIter
%             % Noisify the image
%             for i=1:fppr.nStars
%                 idx = 1 + (i-1)*npsf:i*npsf;
%                 tmp = psfModel(:,idx)*10^(-0.4*mag(k))*F0*F2DN;
%                 psfNoise(:,idx)= poissrnd(tmp + skyBg) - skyBg + eCCD*randn(npsf);
%             end
%             % Cn2h retrieval
%             fppr.psfSky     = psfNoise;
%             fppr.weightInit = ones(1,nL)/nL*atm.r0^(-5/3);
%             fppr            = fppr.getParameters();
%             Wn(:,k,t)      = fppr.xfinal;
%             fW(:,k,t)      = fppr.dweights;
%         end
%     end
%     t0=toc;
%     dWn = 100*bsxfun(@rdivide,abs(bsxfun(@minus,Wn,w0')),w0');
%     
%     fitswrite([dWn,fW],[savePath,'fits/cn2hErrorvMag_',int2str(nL),'layers_ps',int2str(psInMas),'mas.fits']);
% else
%     out = fitsread([savePath,'fits/cn2hErrorvMag_',int2str(nL),'layers_ps',int2str(psInMas),'mas.fits']);
%     dWn = out(:,1:end/2,:);
%     fW  = out(:,end/2+1:end,:);
%     [nL,nMag,nIter] = size(dWn);
%     mag = 0:1:nMag-1;
% end
% 
% % Plot Cn2h error v magnitude
% dW = mean(dWn,3);
% eW = 3*sqrt(var(dWn,[],3)/nIter + median(fW.^2,3));
% %Low altitude bins <1km
% dW1= sqrt(mean(dW(idLow,:).^2,1));
% eW1= sqrt(mean(eW(idLow,:).^2,1));
% %Mid altitude bins <8km
% dW2= sqrt(mean(dW(idMid,:).^2,1));
% eW2= sqrt(mean(eW(idMid,:).^2,1));
% %High altitude bins <16km
% dW3= sqrt(mean(dW(idHigh,:).^2,1));
% eW3= sqrt(mean(eW(idHigh,:).^2,1));
% %High altitude bins <16km
% dW4= sqrt(mean(dW(idVhigh,:).^2,1));
% eW4= sqrt(mean(eW(idVhigh,:).^2,1));
% 
% h = figure;
% plot(mag,dW1,'ks--','MarkerFaceColor','k','MarkerSize',7,'LineWidth',1.2);hold on;
% plot(mag,dW2,'bo--','MarkerFaceColor','b','MarkerSize',7,'LineWidth',1.2);
% plot(mag,dW3,'rd--','MarkerFaceColor','r','MarkerSize',7,'LineWidth',1.2);
% plot(mag,dW4,'mp--','MarkerFaceColor','m','MarkerSize',7,'LineWidth',1.2);
% xlabel('Calibration star magnitude (mag)','interpreter','latex','FontSize',24);
% ylabel('$C_n^2(h)$ error (\%)','interpreter','latex','FontSize',24);
% legend({'$h< 1$km','8km$> h\geq$1km','16km$>h\geq$8km','$h\geq 16$km'},'interpreter','latex',...
%     'FontSize',24,'Location','southeast','AutoUpdate','off');
% psfrTools.makeAxisSquared(h);
% inBetween = [(dW1-eW1) fliplr(dW1+eW1)];
% fill([mag fliplr(mag)], inBetween, 'k','FaceAlpha',.15,'LineStyle','none');
% inBetween = [(dW2-eW2) fliplr(dW2+eW2)];
% fill([mag fliplr(mag)], inBetween, 'b','FaceAlpha',.15,'LineStyle','none');
% inBetween = [(dW3-eW3) fliplr(dW3+eW3)];
% fill([mag fliplr(mag)], inBetween, 'r','FaceAlpha',.15,'LineStyle','none');
% inBetween = [(dW4-eW4) fliplr(dW4+eW4)];
% fill([mag fliplr(mag)], inBetween, 'm','FaceAlpha',.15,'LineStyle','none');
% xlim([mag(1)-1,mag(end)+1]);
% plot([xlim()],[1,1],'k--');
% text(1,1.2,'1\%','interpreter','latex','FontSize',24);
% plot([xlim()],[10,10],'k--');
% text(1,11.5,'10\%','interpreter','latex','FontSize',24);
% set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
% set(gca,'YScale','log')



% close all;
% %7-layers fit
% msk7 = [1292,568,295,111,71,39,19];
% msk6 = [1292,568,295,111,71,29];
% msk5 = [1292,568,295,91,29];
% msk4 = [1292,568,295,95];
% msk3 = [1292,568,138];
% 
% msk = {msk7,msk6,msk5,msk4,msk3};
% nM = length(msk);
% W = zeros(nL,nM);
% 
% for t=1:nM
%     fpp = focalPlaneProfiler(tel,src(msk{t}),ngs(1),atm,otf0,psInMas,otf0,...
%         [],[],'nLayer',nL,'weightInit',w0,'fitWeights',true,'heightInit',h0,...
%         'fitHeight',false,'nptZonal',dk,'idxPhase',validPhase,'flagAnisokinetism',true);
%     
%     psfModel       = fpp.psfModel(w0,fpp.otfKernel);
%     fpp.psfSky     = psfModel;
%     fpp.weightInit = ones(1,nL)/nL*atm.r0^(-5/3);
%     fpp            = fpp.getParameters();
%     W(:,t)         = fpp.xfinal;
% end
% dWt = 100*bsxfun(@rdivide,abs(bsxfun(@minus,W,w0')),w0');
% 
% 
% h = figure;
% semilogy(fliplr(3:7),mean(dWt,1),'bs--','MarkerFaceColor','b','MarkerSize',7);hold on
% xlabel('\# PSF','interpreter','latex','FontSize',24);
% ylabel('$C_n^2(h)$ error (\%)','interpreter','latex','FontSize',24);
% set(gca,'FontSize',24,'FontName','cmr12','TickLabelInterpreter','latex');
% plot([xlim()],[1,1],'k--');
% text(3,1.2,'1\%','interpreter','latex','FontSize',24);
% psfrTools.makeAxisSquared(h)