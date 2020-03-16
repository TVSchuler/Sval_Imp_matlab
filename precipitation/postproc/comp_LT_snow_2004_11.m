%% PROGRAM 
%% WRITTEN TVS, mar 2005

% ----------------------------------------------
% PROGRAM
% ----------------------------------------------
clear all
close all

if(ispc)
    home        = 'M:/';
    datadisk    = 'K:/gggstaff/thomasc/';
%     datadisk    = 'K:/';
else
    home        = '/uio/kant/geo-natg-u1/thomasc/';
    datadisk    = [home,'datadisk/'];
end
% enable asciigridread etc...
addpath(strcat(home,'matlabtools/'));          %% 
% enable SNCTOOLS & mexnc
addpath ([home,'matlabtools/snctools2011/mexcdf/mexnc']);
addpath ([home,'matlabtools/snctools2011/mexcdf/snctools']);
javaaddpath([home,'matlabtools/snctools2011/netcdfAll-4.2.jar']);


%% define filenames
demname2     = [datadisk,'austfonna/dem_sval1km/sval1km_wgs84.txt'];
demname     = [datadisk,'austfonna/dem_sval1km/na_1km_wgs84.txt'];

%%
% years      = [1999 2004 2005 2006 2007];
years      = [2004 2005 2006 2007 2008 2010 2011];
ny        = numel(years);
dens       = [.375,.345,.395,.4,.35,.35,.35];

%% define glacier file
glacname  = strcat(datadisk,'spice/austfonna_dem/austfonna_outline.txt');
glacier   = load(glacname);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ DATASETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read DEM
% import DEM
[zzz, ncol, nrow, xll, yll, cellsize, nodata] = asciigridread(demname);

xkoo        = xll+(0:ncol-1)*cellsize;
ykoo        = yll+(0:nrow-1)*cellsize;
% flip y around: 1rst pnt is upper left corner
ykoo = flipud(ykoo');


% import DEM2
[zz2, ncol2, nrow2, xll2, yll2, cellsize2, nodata2] = asciigridread(demname2);

xko2        = xll2+(0:ncol2-1)*cellsize2;
yko2        = yll2+(0:nrow2-1)*cellsize2;
% flip y around: 1rst pnt is upper left corner
yko2 = flipud(yko2');

[X2,Y2]       = meshgrid(xko2,yko2);
%% load LT_model results

% datafile     = 'results_1406_0411_par_constNm';
% datafile     = 'results_1906_20042011';
% datafile     = 'results_transtau_test';
disp('reading the result grids...')
datafile     = '../results_monthly/LT2014_20150831T215840.mat';
% % datafile     = '../results/results_0307_19792013_tauT_qweight_deg2480'; % results are on the X2,Y2 grid
% % datafile     = 'results_2806_20032013_tauT_qweight_UV3_RH3_deg2480'; % results are on the X2,Y2 grid
% % datafile     = 'results_0107_20032013_tauT_qweight_deg2480'; % results are on the X2,Y2 grid
% % datafile     = 'results_0207_20032013_tauT_qweight_RH80';
% % datafile = 'results_1001_20042011_tauJarosch_urev';
load(datafile);            %%'monthdate','monthprec','X','Y'
disp('done')

[X,Y]       = meshgrid(xkoo,ykoo);% read DEM
% number of years
nY = 8; % 2003/04 2004/05 2005/06 2006/07 2007/08 2009/10 2010/11, 7 years spanning a period of 8 yrs

% make winter prec
dv            = datevec(monthdate);
nyy            = numel(unique(dv(:,1)))-1;%-2; % dataset until mar 2013, so last season not complete

year          = zeros(nyy,1);
winprec       = zeros(nyy,nrow,ncol);
% small hack... because we have overwritten the X,Y grid for monthprec
% (which is flipped)
% yko3 = flipud(yko2);
% [~,ir,ir2]    = intersect(ykoo,yko3);
[~,ir,ir2]    = intersect(ykoo,yko2);
% 
[~,ic,ic2]    = intersect(xkoo,xko2);
% index to accumulation period (Sep-Apr)
idt           = find(dv(:,2)==9);

%%
for i=1:numel(year)
    
    % for some reason, the monthprec matrix is flipud...
    winprec(i,:,:) = sum(monthprec(idt(i):idt(i)+7,ir2,ic2),1)*1e3; % conversion m --> mm
%     winprec(i,:,:) = sum(monthprec(idt(i):idt(i)+7,:,:),1)*1e3; % conversion m --> mm
%     winprec(i,:,:) = sum(monthprec((i-1)*12+(1:8),:,:),1)*1e3; % conversion m --> mm
%     if(i==4) % 2007
%         winprec(i,:,:) = sum(monthprec((i-1)*12+(0:8),ir2,ic2),1)*1e3; % conversion m --> mm
%     end
    % turn winprec around to match orientation of DEM (zzz)
     winprec(i,:,:) = flipud(squeeze(winprec(i,:,:)));
    tmp            = datevec(monthdate(idt(i)+7));
    disp(datestr(tmp))
    year(i)        = tmp(1);
end

% load colormap
load('RE_colors2');  % cmap


snowpath   = '../../../snow/';

% for histograms:
xout   = -500:50:500;
nxout  = length(xout);

n      = zeros(ny,nxout);

for iy     = 1:ny
    disp(num2str(years(iy)))
    idy        = find(year==years(iy));
    if(~isempty(idy))
    tmp        = (squeeze(winprec(idy,:,:)));
    datafile   = strcat(snowpath,'new_snow_',num2str(years(iy)),'_res_1000')
    d{iy}      = load(datafile);
%     m{iy}      = win_TP_19992007(idy,d{iy}.nval);
    m{iy}      = interp2(X,Y,tmp,d{iy}.x_res*1e3,d{iy}.y_res*1e3);
    s{iy}      = d{iy}.s_mean.*dens(iy).*1e3;   % m --> mm
    axmax      = ceil(max([max(s{iy}), max(m{iy})])/500)*500;
    r          = corrcoef(s{iy},m{iy});
    r2 (iy)    = r(1,2).^2;
    bias(iy)   = mean(s{iy}-m{iy});
    rms(iy)    = sqrt(mean((s{iy}-m{iy}).^2));
    stdev(iy)  = std(s{iy}-m{iy});
    
%     if(iy==1)
%         s_v    = s{iy};
%         m_v    = m{iy};
%     else
%         s_v    = [s_v;s{iy}];
%         m_v    = [m_v;m{iy}];
%     end
    
    
axmax=1500;
figure(1)
subplot(2,ceil(ny/2),iy)
plot(s{iy},m{iy},'.')
hold on
plot([0 axmax],[0 axmax],'k')
% axis([0 axmax 0 axmax]);
axis equal
axis square
title(num2str(years(iy)))
xlabel('measured')
ylabel('modeled')
set(gcf,'Position',[300 100 1500 300]);
axis([0 axmax 0 axmax])
text(1.2,0.3,['r\^2 = ',num2str(r2(iy),2)],'units','normalized')
text(1.2,0.2,['bias = ',num2str(bias(iy),5)],'units','normalized')
text(1.2,0.1,['\sigma = ',num2str(stdev(iy),5)],'units','normalized')
 


figure(2)
subplot(2,ceil(ny/2),iy)
n(iy,:) = hist(m{iy}-s{iy},xout);
bar(xout,n(iy,:)/length(m{iy}))
hold on
plot([0 0],[0 0.5],'k')
title(num2str(years(iy)))
xlabel('mod-obs')
ylabel('modeled')
set(gcf,'Position',[300 100 1500 300]);
axis([-300 300 0 0.5])
set(gca,'XTick',-200:100:200);

figure
colormap(jet)
imagesc(xkoo/1e3,ykoo/1e3,squeeze(winprec(idy,:,:)))
hold on
contour(X./1e3,Y./1e3,zzz,0:100:1000,'k')
axis xy equal tight
axis([610 740 8800 8930])
plot(glacier(:,1)/1e3,glacier(:,2)/1e3,'r','linewidth',2)
% plot(d{iy}.x_res,d{iy}.y_res,'.k')
scatter(d{iy}.x_res,d{iy}.y_res,25,s{iy},'filled','markeredgecolor','k')
title(num2str(years(iy)))
caxis([0 1500])
colorbar

figure
% imagesc(xkoo/1e3,ykoo/1e3,squeeze(winprec(idy,:,:)))
contour(X./1e3,Y./1e3,zzz,0:100:1000,'k')
axis xy equal tight
hold on
plot(glacier(:,1)/1e3,glacier(:,2)/1e3,'r','linewidth',2)
axis([610 740 8800 8930])
% plot(d{iy}.x_res,d{iy}.y_res,'.k')
scatter(d{iy}.x_res,d{iy}.y_res,25,(m{iy}-s{iy})./s{iy},'filled')%,'markeredgecolor','k')
title(num2str(years(iy)))
caxis([-.5 .5])
colormap(cmap)
colorbar

figure(18)
plot(s{iy},m{iy},'.')
hold all
axis([0 1500 0 1500])
axis square

end

end

% title(exp_nr)
text(0.7,0.3,['r\^2 = ',num2str(mean(r2),2)],'units','normalized')
text(0.7,0.25,['bias = ',num2str(mean(bias),5)],'units','normalized')
text(0.7,0.2,['\sigma = ',num2str(mean(stdev),5)],'units','normalized')
legend('2004','2005','2006','2007','2008','2010','2011')
plot([0 1500],[0 1500],'k')



% r_t        = corrcoef(s_v,m_v);
% r_t2       = r(1,2).^2;
% bias_t     = mean(s_v-m_v);
% rms_t      = sqrt(mean((s_v-m_v).^2));

% figure
% plot(s_v,m_v,'.')
% axis([0 1500 0 1500])
% hold on
% plot([0 1500],[0 1500],'k')
% axis square
% % title(exp_nr)
% text(0.7,0.3,['rsq = ',num2str(r_t2,2)],'units','normalized')
% text(0.7,0.25,['bias = ',num2str(bias_t,5)],'units','normalized')
% text(0.7,0.2,['rms = ',num2str(rms_t,5)],'units','normalized')
% 

% weighted histograms...
for iy = 1:ny
    nb(iy,:)=n(iy,:)/length(m{iy});
end

nbb = mean(nb,1);


figure
% [na,xouta] = hist(m_v-s_v,-250:50:500);
% bar(xouta,na)
bar(xout,nbb)
hold on
plot([0 0],[0 0.5],'k')
xlim([-300 300])
ylim([0 0.3])
%% read parameter values from runloop_ecmwf_newthermo.scr
% exp_par=importdata([nc_path,'/runloop_ecmwf_newthermo.scr'],'=',44);


% %% result structure for output
% results    = struct('exp_ID',{exp_nr},'Nm',{exp_par.data(1)},'tau',{exp_par.data(2)},'years',{years},'rsq',{r2},...
%     'bias',{bias},'rms',{rms},'rsq_all',{r_t2},'bias_all',{bias_t},'rms_all',{rms_t});
% 
% outfile    = ['./results_snow/',exp_nr,'_snow_results']
% % save(outfile,'results');