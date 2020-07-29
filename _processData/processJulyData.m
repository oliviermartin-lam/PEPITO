%% MANAGE workspaces
close all; % close all figures
clear all; % clear the workspace
clc;

path_oomao = '/home/omartin/Projects/SIMULATIONS/OOMAO/mfiles/'; % the oomao path
path_workspace = '/home/omartin/Projects/PEPITO/'; % the simulation folder path
path_res= '/home/omartin/Projects/PEPITO/Results/'; % the simulation folder path

addpath(genpath(path_oomao),genpath(path_workspace));

%% GET DATA ID
path_data = '/run/media/omartin/OlivierMartinHDD/DATA/PEPITO_DATA/'; % the simulation folder path
folder_data = dir(path_data);
dataID = {folder_data(contains({folder_data.name},'.fits')).name};
nObj = numel(dataID) - 6;

%% INSTANTIATING PEPITO
fitCn2 = true;
pep = pepito('parFile_pt5m','fitCn2',fitCn2);

%% LOOP ON DATA
if fitCn2
    res = zeros(2,17,nObj);
else
    res = zeros(2,10,nObj);
end
method = 'PSF';

for k=1:1
    fprintf('I am processing the %dth file\n',k);
    % READ DATA
    path_im = [path_data,cell2mat(dataID(k+6))];
    hdr = fitsinfo(path_im);hdr = hdr.PrimaryData.Keywords;
    
    % GET ACQUISITION TIME
    h = cell2mat( hdr(strcmp(hdr(:,1),'UTHOUR') , 2) );
    m = cell2mat( hdr(strcmp(hdr(:,1),'UTMIN') , 2) );
    s = cell2mat( hdr(strcmp(hdr(:,1),'UTSEC') , 2) );
    uttime = h + m/60 + s/3600;
    utday = cell2mat( hdr(strcmp(hdr(:,1),'UTDAY') , 2) );
    utmonth = cell2mat( hdr(strcmp(hdr(:,1),'UTMONTH') , 2) );
    utyear = cell2mat( hdr(strcmp(hdr(:,1),'UTYEAR') , 2) );
    
    % GET EXPOSURE
    tExp = cell2mat( hdr(strcmp(hdr(:,1),'EXP') , 2) );
    conjAlt = cell2mat( hdr(strcmp(hdr(:,1),'CONJALT') , 2) );
    
    imSE = fitsread(path_im);
    if ndims(imSE) == 3
        pep.runPepito(imSE,'fitL0',false,'method',method);
        
        % SAVING RESULTS
        %res(1,:,k) = [utyear,utmonth,utday,uttime,tExp,conjAlt,pep.r0,pep.Cn2h,pep.theta,pep.vib];
        %res(2,:,k) = [zeros(1,6),pep.dr0,pep.dCn2h,pep.dtheta,reshape(pep.dvib,[],3)];
    end    
end
%fitswrite(res,[path_res,'pepitoCn2Estimation_July_',method,'.fits']);

%% display
close all;
if strcmp(method,'PSF')
    pep.displayResults()
else
    x = linspace(0,(pep.nBox-1)/2,pep.nBox/2) * pep.psInMas/1e3;
    fontsize = 20;
    clf;
    plot(x,pep.ydata,'bs','MarkerFaceColor','b','MarkerSize',7);
    hold on;
    plot(x,pep.mod(pep.xfinal,pep.xdata),'k--','linewidth',1.5);
    xlabel('Angular separation (arcsec)','interpreter','latex','fontsize',fontsize);
    ylabel('PSF profile','interpreter','latex','fontsize',fontsize);    
    legend({'On-sky PSF','Fitted Model'},'interpreter','latex','fontsize',fontsize);
    set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
    xlim([0,max(x)])
    pbaspect([1,1,1]);
end


%% Display
close all;
fontsize = 24;
% Grab results
res = fitsread([path_res,'pepitoCn2Estimation_July_PSF.fits']);

r0 = 1e2*squeeze(res(1,7,:));
dr0 = 1e2*squeeze(res(2,7,:));
th0 = squeeze(res(1,14,:));
dth0 = squeeze(res(2,14,:))/1e2;
vib = squeeze(res(1,8:end,:));
dvib = squeeze(res(2,8:end,:));
utyear = squeeze(res(1,1,:));utmonth = squeeze(res(1,2,:));
utday = squeeze(res(1,3,:));uttime = squeeze(res(1,4,:));


dd = unique(utday);
nDay = numel(dd);

id = find(r0 > 4);
figure;
errorbar(uttime(id),r0(id),dr0(id),'bs','MarkerFaceColor','b','MarkerSize',7);    
xlabel('UT time (h)','interpreter','latex','fontsize',fontsize)
ylabel('Line of sight $r_0$ at 500 nm (cm)','interpreter','latex','fontsize',fontsize)
set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
pbaspect([1,1,1]);

id = find(th0 < 2.5);
figure;
errorbar(uttime(id),th0(id),dth0(id),'bs','MarkerFaceColor','b','MarkerSize',7);    
xlabel('UT time (h)','interpreter','latex','fontsize',fontsize)
ylabel('Line of sight $\theta_0$ at 500 nm (cm)','interpreter','latex','fontsize',fontsize)
set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
pbaspect([1,1,1]);

figure;
plot(uttime,sqrt(vib(1,:)),'bs','MarkerFaceColor','b','MarkerSize',7);    
hold on;
plot(uttime,sqrt(vib(2,:)),'ro','MarkerFaceColor','r','MarkerSize',7);    
xlabel('UT time (h)','interpreter','latex','fontsize',fontsize)
ylabel('Jitter (arcsec)','interpreter','latex','fontsize',fontsize)
legend({'X-axis','Y-axis'},'interpreter','latex','fontsize',fontsize);
set(gca,'FontSize',fontsize,'FontName','cmr12','TickLabelInterpreter','latex');
pbaspect([1,1,1]);