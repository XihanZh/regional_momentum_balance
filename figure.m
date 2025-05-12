%--- figure 1 model domain 

load chp3new.mat
load ekenew.mat
% B2=[230,460,125,343];%large meander
B2new=[230,460,70,343];%large meander


dzf0=repmat(dzf,[1,600,400]);
dzf3d=permute(dzf0,[2,3,1]);
dzf3w=dzf3d.*hacw;
udz=squeeze(sum(ubar.*dzf3w,3,'omitnan'));% 600*400
udzdy=udz.*dyg;
psibt=cumsum(udzdy,2)/1e6;
% d=rdmds('Depth');

x1=151;
x2=401;
tc=20:10:150;
[X,Y]=meshgrid(x,y);


close all
x0=10;
y0=10;
width=800;
height=2000;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(2,1,1);
imagesc(x,y,depth'),axis xy
colormap(ax(1),"parula")
colorbar
clim([0 6000])
clb = colorbar('Ticks',0:1000:6000);
clb.FontSize = 12;
clb.Label.String = 'Bathymetry [m]';
clb.Label.FontSize = 12;
tobj=title('(a)','fontsize',14);
tobj.Position=[x(1)+0.4,y(end)+0.1];
% set(gca,'fontsize',14)
% ylabel(colorbar,'Bathymetry [m]','FontSize',14)
set(ax(1),'XTick',[],'YTick',[])
h1=axes('position',get(ax(1),'position'),'color','none','fontsize',14);
hold on
contour(X,Y,psibt',tc,'LineWidth',2,'ShowText','on')
colormap(h1,[0,0,0])
hold on
plot(x(x1)*ones(size(y(2:end-1))),y(2:end-1),'r--','LineWidth',2)
plot(x(x2)*ones(size(y(2:end-1))),y(2:end-1),'r--','LineWidth',2)
text(x(x1),y(20),'A','fontsize',24)
text(x(x2),y(20),'B','fontsize',24)

% text(x(191),y(50),'Southeast Indian Ridge','color','k','fontsize',14)
% text(x(477),y(340),'Macquarie Ridge','color','k','fontsize',14)

text(x(191),y(50),'Southeast Indian Ridge','color',[1,1,1],'fontsize',16)
text(x(477),y(340),'Macquarie Ridge','color',[1,1,1],'fontsize',16)


ax(2)=subplot(2,1,2);
% imagesc(x,y,ekemap'),axis xy
imagesc(x,y,squeeze(eke(:,:,1))'),axis xy
colorbar
set(gca,'fontsize',14)
clim([0 0.2])
colormap(ax(2),'hot')
set(ax(2),'xtick',[],'ytick',[]);
ylabel(colorbar,'EKE [m^{2}s^{-2}]')


[X,Y]=meshgrid(x,y);
h2=axes('position',get(ax(2),'position'),'color','none','fontsize',14);
hold on
[cn,hn]=contour(X,Y,psibt',tc,'showtext','on','LineWidth',2,'parent',h2);
clabel(cn,hn,0:10:160,'color','c','fontsize',10)
colormap(h2,[1,1,1])
rr=B2new;
plot(x(rr(1))*ones(1,length(rr(3):rr(4))),y(rr(3):rr(4)),'b','LineWidth',2)
plot(x(rr(2))*ones(1,length(rr(3):rr(4))),y(rr(3):rr(4)),'b','LineWidth',2)
plot(x(rr(1):rr(2)),y(rr(3))*ones(1,length(rr(1):rr(2))),'b','LineWidth',2)
plot(x(rr(1):rr(2)),y(rr(4))*ones(1,length(rr(1):rr(2))),'b','LineWidth',2)
tobj=title('(b)','fontsize',14);
tobj.Position=[x(1)+0.4,y(end)+0.1];
set(gca,'fontsize',14)
print -dpng modeldomain.png

%-- figure 2 momentum profile 
load('EPrealheatflux1.mat')
rac3d = repmat(rac,[1,1,200]);
total_area = sum(rac,'all');
water_area = sum(squeeze(sum(rac3d.*hacc,1)),1);
rho0 = 999.8;
load mmtprofile.mat
z1=16; % 80m
z2 = 98; % 2km
geocoriolis = -dphipf;
ageo = corispf-geocoriolis;

% Create figure and set size
x0 = 10;
y0 = 10;
width = 800;
height = 3200;
set(gcf, 'Position', [x0, y0, width, height]);

% Define normalized positions (between 0 and 1)
topPos = [0.1, 0.55, 0.7, 0.4];   % [left, bottom, width, height]
botPos = [0.1, 0.05, 0.7, 0.4];
sidePos = [0.81, 0.05, 0.15, 0.4]; % Narrow panel beside lower one


ax(1) = subplot('Position', topPos);
plot(rho0*advpf(1:z1),zc(1:z1),'LineWidth',2)
hold on
plot(rho0*ageo(1:z1),zc(1:z1),'LineWidth',2)
plot(rho0*geocoriolis(1:z1),zc(1:z1),'LineWidth',2)
plot(rho0*disspf(1:z1),zc(1:z1),'LineWidth',2)
plot(rho0*dphipf(1:z1),zc(1:z1),'LineWidth',2)
plot(rho0*vsflxtotpf(1:z1)+rho0*extpf(1:z1),zc(1:z1),'c','LineWidth',2)
grid on
xlim([-6e9 6e9])
ylim([-80 0])
ylabel('Depth [m]')
set(ax(1),'fontsize',14)
title('(a)','position',[-5.8e9 0.1])



ax(2) = subplot('Position', botPos);
plot(rho0*advpf,zc,'LineWidth',2)
hold on

plot(rho0*ageo,zc,'LineWidth',2)
plot(rho0*geocoriolis,zc,'LineWidth',2)
plot(rho0*disspf,zc,'LineWidth',2)
plot(rho0*dphipf,zc,'LineWidth',2)
plot(rho0*vsflxtotpf+rho0*extpf,zc,'c','LineWidth',2)
xlim([-6e9 6e9])
ylim([zc(end) 0])
grid on
title('(b)','position',[-5.8e9 0.1])
ylabel('Depth [m]')
xlabel('zonal momentum terms [N m^{-1}]')
set(ax(2),'fontsize',14,'Layer','top')
plot(linspace(-6*1e9,6*1e9,20),ones(1,20)*zc(z1),'k--')
plot(linspace(-6*1e9,6*1e9,20),ones(1,20)*zc(z2),'m--')
legend({'advection','ageostrophic Coriolis','geostrophic Coriolis','horizontal & bottom friction',...
    'PGF','vertical friction','80m','2km'},'fontsize',14,'Location','best')


ax(3) = axes('Position', sidePos);
plot((total_area-water_area)/total_area,zc,'LineWidth',2)
grid on
set(ax(3),'fontsize',14,'Layer','top')
title('(c)','position',[0.1 0.1])
yticks('')

% ---- figure 3 momentum maps interior 
load mmtmapsrealcomp.mat
load chp3.mat
load blue_red_saturated.mat
D=[46,600-45,46,400-45];
rr = D;
rho0 = 999.8;

close all
x0=10;
y0=10;
width=800;
height=4000;
set(gcf,'position',[x0,y0,width,height])
% fig=gcf;

%-- middle depth
ax(1)=subplot(4,1,1);
imagesc(x(rr(1):rr(2)),y(rr(3):rr(4)),rho0*coris_mapmid')
colormap(map)
% colorbar
clim(ax(1),[-25 25])
% colorbar('Ticks',-0.02:0.01:0.02)
cl1 = colorbar('Ticks',-20:10:20);
tobj=title('(a)','fontsize',14);
tobj.Position=[x(46)+0.4,y(40)];
text(x(280),y(20),'Coriolis','fontsize',14)
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
% ylabel(colorbar,'[N/m^{2}]')
ylabel(cl1,'[N/m^{2}]')
set(gca,'fontsize',14)


ax(2)=subplot(4,1,2);
imagesc(x(rr(1):rr(2)),y(rr(3):rr(4)),rho0*dphi_mapmid')
colormap(map)
colorbar
clim(ax(2),[-25 25])
% colorbar('Ticks',-0.02:0.01:0.02)
cl2 = colorbar('Ticks',-20:10:20);
% title('PGF interior')
tobj=title('(b)','fontsize',14);
tobj.Position=[x(46)+0.4,y(40)];
text(x(280),y(20),'PGF','fontsize',14)

yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
% ylabel(colorbar,'[N/m^{2}]')
ylabel(cl2,'[N/m^{2}]')
set(gca,'fontsize',14)

ax(3)=subplot(4,1,3);
imagesc(x(rr(1):rr(2)),y(rr(3):rr(4)),rho0*adv_mapmid')
colormap(map)
colorbar
clim(ax(3),[-5e-1,5e-1])
cl3 = colorbar('Ticks',-5e-1:2.5e-1:5e-1);
tobj=title('(c)','fontsize',14);
tobj.Position=[x(46)+0.4,y(40)];
text(x(280),y(20),'advection','fontsize',14)
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel(cl3,'[N/m^{2}]')
set(gca,'fontsize',14)


ax(4)=subplot(4,1,4);
imagesc(x(rr(1):rr(2)),y(rr(3):rr(4)),rho0*diss_mapmid')
colormap(map)
% colorbar
clim(ax(4),[-5e-2,5e-2])
cl4 = colorbar('Ticks',-5e-2:2.5e-2:5e-2);
tobj=title('(d)','fontsize',14);
tobj.Position=[x(46)+0.4,y(40)];
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel(cl4,'[N/m^{2}]')
text(x(280),y(20),'friction','fontsize',14)
set(gca,'fontsize',14)
print -dpng mmtmaps_midnew.png

% --- figure 4 momentum maps bottom layer
close all
x0=10;
y0=10;
width=800;
height=4000;
set(gcf,'position',[x0,y0,width,height])
% fig=gcf;


ax(1)=subplot(4,1,1);
imagesc(x(rr(1):rr(2)),y(rr(3):rr(4)),rho0*coris_mapbot')
colormap(map)
% colorbar
clim(ax(1),[-25 25])
cl1 = colorbar('Ticks',-20:10:20);
tobj=title('(a)','fontsize',14);
tobj.Position=[x(46)+0.4,y(40)];
text(x(280),y(20),'Coriolis','fontsize',14)
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel(cl1,'[N/m^{2}]')
set(gca,'fontsize',14)


ax(2)=subplot(4,1,2);
imagesc(x(rr(1):rr(2)),y(rr(3):rr(4)),rho0*dphi_mapbot')
colormap(map)
clim(ax(2),[-25 25])
cl2 = colorbar('Ticks',-20:10:20);
tobj=title('(b)','fontsize',14);
tobj.Position=[x(46)+0.4,y(40)];
text(x(280),y(20),'PGF','fontsize',14)
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel(cl2,'[N/m^{2}]')
set(gca,'fontsize',14)

ax(3)=subplot(4,1,3);
imagesc(x(rr(1):rr(2)),y(rr(3):rr(4)),rho0*adv_mapbot')
colormap(map)
clim(ax(3),[-5e-1,5e-1])
cl3 = colorbar('Ticks',-5e-1:2.5e-1:5e-1);
% title('friction interior')
tobj=title('(c)','fontsize',14);
tobj.Position=[x(46)+0.4,y(40)];
text(x(280),y(20),'advection','fontsize',14)
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel(cl3,'[N/m^{2}]')
set(gca,'fontsize',14)



ax(4)=subplot(4,1,4);
imagesc(x(rr(1):rr(2)),y(rr(3):rr(4)),rho0*diss_mapbot')
colormap(map)

clim(ax(4),[-5e-1,5e-1])
cl4 = colorbar('Ticks',-5e-1:2.5e-1:5e-1);
% title('friction bottom')
tobj=title('(d)','fontsize',14);
tobj.Position=[x(46)+0.4,y(40)];
text(x(280),y(20),'friction','fontsize',14)
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel(cl4,'[N/m^{2}]')
set(gca,'fontsize',14)
print -dpng mmtmaps_botnew.png


%--- figure 5 TFS maps 
xhr = load('TFShighres.mat','x').x;
yhr = load('TFShighres.mat','y').y;
ubar = load('EPrealheatflux1.mat','ubar').ubar;
psibthr = load('psibthighres.mat').psibt;
load TFShighres.mat
load TFSnewmode.mat
load blue_red_saturated.mat
load chp3new.mat

dzf0=repmat(dzf,[1,600,400]);
dzf3d=permute(dzf0,[2,3,1]);
dzf3w=dzf3d.*hacw;
udz=squeeze(sum(ubar.*dzf3w,3,'omitnan'));% 600*400
udzdy=udz.*dyg;
psibt=cumsum(udzdy,2)/1e6;
tct=0:10:160;
[X,Y]=meshgrid(x,y);

% B2=[230,460,125,343];%large meander
B2new=[230,460,70,343];
x0=10;
y0=10;
width=800;
height=2000;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(2,1,1);
imagesc(x(2:end-1),y(2:end),TFSreg'),axis xy
colormap(ax(1),map1)
clim([-2 2])
colorbar('Ticks',-2:0.5:2,'fontsize',12)
tobj=title('(a)','fontsize',14);
tobj.Position=[x(1)+0.4,y(end)+0.1];
% set(gca,'fontsize',14)
ylabel(colorbar,'TFS [N/m^{2}]','FontSize',12)
set(ax(1),'XTick',[],'YTick',[])
h1=axes('position',get(ax(1),'position'),'color','none','fontsize',14);
hold on
contour(X,Y,psibt',tct,'LineWidth',2,'ShowText','on')
colormap(h1,[0,0,0])
hold on
rr=B2new;
plot(x(rr(1))*ones(1,length(rr(3):rr(4))),y(rr(3):rr(4)),'r','LineWidth',2)
plot(x(rr(2))*ones(1,length(rr(3):rr(4))),y(rr(3):rr(4)),'r','LineWidth',2)
plot(x(rr(1):rr(2)),y(rr(3))*ones(1,length(rr(1):rr(2))),'r','LineWidth',2)
plot(x(rr(1):rr(2)),y(rr(4))*ones(1,length(rr(1):rr(2))),'r','LineWidth',2)
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})

Bhigh = [900,1800,400,1200];
rr=Bhigh;
ax(2)=subplot(2,1,2);
imagesc(xhr(2:end-1),yhr(2:end-1),TFShighres'),axis xy
colorbar('Ticks',-2:0.5:2,'fontsize',12)
set(ax(2),'xtick',[],'ytick',[]);
% set(gca,'fontsize',14)
clim([-2 2])
colormap(ax(2),map1)
% set(ax(2),'xtick',[],'ytick',[]);
ylabel(colorbar,'TFS [N/m^{2}]','FontSize',12)

[Xh,Yh]=meshgrid(xhr(2:end-1),yhr(2:end-1));
h2=axes('position',get(ax(2),'position'),'color','none','fontsize',14);
hold on
% [cn,hn]=contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),etabar(2:end-1,2:end-1)',-0.5:0.1:0.5,'showtext','on','LineWidth',2,'parent',h2);
% clabel(cn,hn,-0.5:0.1:0.5,'color','c','fontsize',10)
[cn,hn]=contour(Xh,Yh,psibthr(2:end-1,2:end-1)',0:10:160,'showtext','on','LineWidth',2,'parent',h2);
clabel(cn,hn,0:10:160,'color','k','fontsize',10)
colormap(h2,[0,0,0])
hold on
plot(xhr(rr(1))*ones(1,length(rr(3):rr(4))),yhr(rr(3):rr(4)),'r','LineWidth',2)
plot(xhr(rr(2))*ones(1,length(rr(3):rr(4))),yhr(rr(3):rr(4)),'r','LineWidth',2)
plot(xhr(rr(1):rr(2)),yhr(rr(3))*ones(1,length(rr(1):rr(2))),'r','LineWidth',2)
plot(xhr(rr(1):rr(2)),yhr(rr(4))*ones(1,length(rr(1):rr(2))),'r','LineWidth',2)
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
tobj=title('(b)','fontsize',14);
tobj.Position=[xhr(1)+0.4,yhr(end)+0.1];

set(gca,'fontsize',14)
print -dpng TFSmap_Whole.png

%--- figure 6 EP flux section 
x1=151;
x2=401;
y = yc(1,:);
x = xc(:,1);
[X1,Z1]=meshgrid(x(2:end),zc);
[X2,Z2]=meshgrid(x(2:end-1),zf(2:end-1));
[X3,Z3]=meshgrid(x(2:end-1),zc(2:end-1));

lat=290;

figure
x0=10;
y0=10;
width=800;
height=3200;
set(gcf,'position',[x0,y0,width,height])

ax(2)=subplot(3,1,1);
pcolor(X1,Z1,log10(squeeze(ekereal(:,lat,:)))')
shading flat
colorbar
% xlabel('Longitude')
ylabel('Depth [m]')
% title('EKE')
% title(strcat('y=',num2str(y(i))),'fontsize',14)
title('(a)','position',[135.5 0.1],'fontsize',14)
colormap(ax(2),jet)
set(gca,'fontsize',14,'Layer','top')
ylabel(colorbar,'EKE [m^{2}s^{-2}]','fontsize',14)
clim([-4.5 -1])
hold on
plot(x(x1)*ones(size(zc)),zc,'k--','LineWidth',2)
plot(x(x2)*ones(size(zc)),zc,'k--','LineWidth',2)
% plot(x(201)*ones(size(zc)),zc,'k--','LineWidth',2)
% plot(x(470)*ones(size(zc)),zc,'k--','LineWidth',2)

ax(3)=subplot(3,1,2);
pcolor(X2,Z2,squeeze(crossep(:,lat,:))')
shading flat
colorbar
clim([-1.5e-3 1.5e-3])
hold on
plot(x(x1)*ones(size(zf(2:end-1))),zf(2:end-1),'k--','LineWidth',2)
plot(x(x2)*ones(size(zf(2:end-1))),zf(2:end-1),'k--','LineWidth',2)
% plot(x(201)*ones(size(zf(2:end-1))),zf(2:end-1),'k--','LineWidth',2)
% plot(x(470)*ones(size(zf(2:end-1))),zf(2:end-1),'k--','LineWidth',2)
% xlabel('Longitude')
ylabel('Depth [m]')
% title('EP flux')
title('(b)','position',[135.5 0.1],'fontsize',14)
ylabel(colorbar,'[m^2 s^{-2}]')
colormap(ax(3),jet)
set(gca,'fontsize',14,'Layer','top')

ax(4)=subplot(3,1,3);
pcolor(X3,Z3,squeeze(crossepdz(:,lat,:))')
shading flat
colorbar
% caxis([-2.5e-4 2.5e-4])
clim([-1e-6 1e-6])
hold on
plot(x(x1)*ones(size(zc(2:end-1))),zc(2:end-1),'k--','LineWidth',2)
plot(x(x2)*ones(size(zc(2:end-1))),zc(2:end-1),'k--','LineWidth',2)
% plot(x(201)*ones(size(zc(2:end-1))),zc(2:end-1),'k--','LineWidth',2)
% plot(x(470)*ones(size(zc(2:end-1))),zc(2:end-1),'k--','LineWidth',2)
colormap(ax(4),jet)
xlabel('Longitude')
ylabel('Depth [m]')
ylabel(colorbar,'[m s^{-2}]')
% title('vertical derivative of EP flux')
title('(c)','position',[135.5 0.1],'fontsize',14)
set(gca,'fontsize',14,'Layer','top')
print -dpng EPfluxsection.png

% figure 7 EP flux profiles
ubar0 = load('saltempchp3newmodel.mat','ubar0').ubar0;

dyg3d=repmat(dyg,[1,1,200]);
msk= psibt>=40 & psibt <=130;
msk3d=repmat(msk,[1,1,200]);
% transport=squeeze(sum(ubar0.*dyg3d.*dzf3w.*msk3d,2));
transport=squeeze(sum(ubar0.*dyg3d.*msk3d,2));
msk0 = psibt(2:end-1,2:end)>40 & psibt(2:end-1,2:end)<130;
msk3d0=repmat(msk0,[1,1,199]);% 599*399*199
rac3d=repmat(rac(2:end-1,2:end),[1,1,199]);% 599*399*199
rac3d(isnan(crossep))=nan; 

rac3dmsk=rac3d.*msk3d0;
crossepdn=crossep(x1:x2,:,:).*rac3dmsk(x1:x2,:,:);
area2=sum(squeeze(sum(rac3dmsk(x1:x2,:,:),1,'omitnan')),1,'omitnan');
area2(area2==0)=nan;
crossepdnavg=sum(squeeze(sum(crossepdn,1,'omitnan')),1,'omitnan')./area2;

crossepdnavg(1:52)=nan;

figure
x0=10;
y0=10;
width=500;
height=500;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(1,2,1);
plot(transport(x1,:),zc,'LineWidth',2)
hold on
plot(transport(x2,:),zc,'LineWidth',2)
grid on
legend({'A','B'},'fontsize',14,'Location','best')
ylim([-4500 0])
ylabel('Depth [m]')
xlabel('[m^{2}s^{-1}]')
title('(a)','position',[2.5e3 0.1])
set(gca,'fontsize',14)

ax(2)=subplot(1,2,2);
plot(crossepdnavg,zc(2:end),'LineWidth',2)
hold on
plot(EP_fluxVdnavg,zc(2:end),'k--','LineWidth',2)
grid on
ylim([-4500 0])
% xlim([0 2e-4])
xlabel('[m^{2}s^{-2}]')
title('(b)','position',[1e-5 0.1])
set(gca,'fontsize',14)
legend({'cross-stream','meridional'},'fontsize',14,'Location','best')
print -dpng EPfluxprofile.png


