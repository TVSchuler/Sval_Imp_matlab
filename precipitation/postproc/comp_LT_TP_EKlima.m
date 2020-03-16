%% PROGRAM 
%% WRITTEN TVS, mar 2005

% ----------------------------------------------
% PROGRAM
% ----------------------------------------------
% clear all
% close all

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ DATASETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read DEM
% import DEM
[zzz, ncol, nrow, xll, yll, cellsize, nodata] = asciigridread(demname2);

xkoo        = xll+(0:ncol-1)*cellsize;
ykoo        = yll+(0:nrow-1)*cellsize;
% flip y around: 1rst pnt is upper left corner
ykoo = flipud(ykoo');


% Y = flipud(Y);

%% load LT_model results

% datafile     = 'results_1406_0411_par_constNm';
% datafile     = 'results_1906_20042011';
% datafile     = 'results_transtau_test';
disp('reading the result grids...e4')
% % datafile     = 'results_2506_19792013_tauT_qweight'; % results are on the X2,Y2 grid
% % datafile     = '../results/results_0107_20032013_tauT_qweight_deg2480'; % results are on the X2,Y2 grid
% datafile     = '../results/results_0307_19792013_tauT_qweight_deg2480'; % results are on the X2,Y2 grid
% datafil2     = '../results/results_0307_19572002_tauT_qweight_deg2480'; % results are on the X2,Y2 grid
% datafil3     = '../results/results_2009_19792013_TPonly'; % results are on the X2,Y2 grid
% datafil4     = '../results/results_2009_19572002_TPonly'; % results are on the X2,Y2 grid

% datafile = 'K:\gggruppe-prosj\glasio-data\svalbard\austfonna\precipLTplain\ei\LTcluster_ei_Svalbard_monthly';
datafile = 'K:\gggruppe-prosj\glasio-data\svalbard\austfonna\precipLT2015\ei\LT2015_ei_Svalbard_monthly';
datafil2 = 'K:\gggruppe-prosj\glasio-data\svalbard\austfonna\precipLT2015\e4\LT2015_e4_Svalbard_monthly';
datafil3 = 'K:\gggruppe-prosj\glasio-data\svalbard\austfonna\precipLTcluster\ei\TP_ei_Svalbard_monthly_threshold';
% datafile = 'results_1001_20042011_tauJarosch_urev';

e4LT = load(datafil2);

% eiLT.monthprec = monthprec; eiLT.monthdate=monthdate;
% eiTP = load(datafil2);            %% monthdate, monthTP, lat, lon
% eiTP2 = load(datafil3);

% e4LT = load(datafil2);
% eiTP = load(datafil3);            %%'monthdate','monthprec'
% e4TP = load(datafil4);
%%
disp('done')
[X,Y]       = meshgrid(xkoo,ykoo);% use coordinates of the DEM, not the X,Y of the results (which are flipped!!!

% [LON,LAT]        = meshgrid(eiLT.lon,eiLT.lat);
% [TPX,TPY]        = deg2utm33(LAT(:),LON(:));
% TPX              = reshape(TPX,size(LON));
% TPY              = reshape(TPY,size(LON));

%% load e-klima dataset
load('Monthly_T_P_Svalbard'); % station: ID, name, lat, lon, x, y, y, z, date, temp, prec

%% interpolate LT/TP to station location
% era-interim period

pre_ei = find(e4LT.monthdate<datenum(1979,1,1));
nt_e4 = length(e4LT.monthdate(pre_ei));

%% station = struct([]);
disp('interpolating...e4')

for ii = 1:numel(station)
    
    for it=1:nt_e4
%         station{ii}.eiTP(it)=griddata(TPX,TPY,squeeze(eiTP.monthTP(it,:,:)),station{ii}.x,station{ii}.y);
%         station{ii}.eiTP2(it)=griddata(TPX,TPY,squeeze(eiTP2.monthTP(it,:,:)),station{ii}.x,station{ii}.y);
        station{ii}.eiLT(it)=interp2(X,Y,squeeze(e4LT.monthprec(it,:,:)),station{ii}.x,station{ii}.y);
    end
    station{ii}.ei_monthdate = floor(e4LT.monthdate(pre_ei));
end
disp('ferdig')
clear e4LT


disp('reading the result grids...e4')
eiLT = load(datafile);            %%'monthdate','monthprec'
disp('done')

nt_ei = size(eiLT.monthprec);
disp('interpolating...ei')
for ii = 1:numel(station)
    
    for it=1:nt_ei(1)
%         station{ii}.eiTP(it)=griddata(TPX,TPY,squeeze(eiTP.monthTP(it,:,:)),station{ii}.x,station{ii}.y);
%         station{ii}.eiTP2(it)=griddata(TPX,TPY,squeeze(eiTP2.monthTP(it,:,:)),station{ii}.x,station{ii}.y);
        station{ii}.eiLT(nt_e4+it)=interp2(X,Y,squeeze(eiLT.monthprec(it,:,:)),station{ii}.x,station{ii}.y);
    end
    station{ii}.ei_monthdate = [station{ii}.ei_monthdate;floor(eiLT.monthdate)];
end
% era-40 period
% nt_e4 = size(e4TP.monthprec);
% for ii = 1:numel(station)
%     for it=1:nt_e4(1)
% %         station{ii}.e4TP(it)=interp2(X,Y,squeeze(e4TP.monthprec(it,:,:)),station{ii}.x,station{ii}.y);
%         station{ii}.e4LT(it)=interp2(X,Y,squeeze(e4LT.monthprec(it,:,:)),station{ii}.x,station{ii}.y);
%     end
%     station{ii}.e4_monthdate = e4TP.monthdate;
% end

%%
for i=2:5
    figure(i)
    stairs(station{i}.date,station{i}.prec)
    hold all
    stairs(station{i}.ei_monthdate,station{i}.eiLT.*1e3)
% %     stairs(station{i}.ei_monthdate,station{i}.eiTP.*1e3)
%     stairs(station{i}.ei_monthdate,station{i}.eiTP2.*1e3)
% %     stairs(station{i}.ei_monthdate,(station{i}.eiLT-station{i}.eiTP+station{i}.eiTP2).*1e3)
    
    datetick
    title(station{i}.name)
    
    [cd,ia,ib]=intersect(floor(station{i}.date),station{i}.ei_monthdate);
    
    figure(10+i)
%     plot(station{i}.prec(ia).*1.5,station{i}.eiTP(ib).*1e3,'.')
%     hold all
%     plot(station{i}.prec(ia).*1.5,(station{i}.eiLT(ib)-station{i}.eiTP(ib)+station{i}.eiTP2(ib)).*1e3,'.')
    plot(station{i}.prec(ia)*1.5,station{i}.eiLT(ib).*1e3,'.')
    hold all
    xl=get(gca,'XLim');yl=get(gca,'YLim');
    ttt=max([xl,yl]);
    plot([0,ttt],[0,ttt],'k')
    title(station{i}.name)
    axis equal square
    bias(i) = nanmean(station{i}.prec(ia)*1.5-station{i}.eiLT(ib)'.*1e3);
    rms(i) = sqrt(nanmean((station{i}.prec(ia)*1.5-station{i}.eiLT(ib)'.*1e3).^2));
    ns(i) = length(cd)
    mobs(i) = nanmean(station{i}.prec(ia));
    dv = datevec(cd);
    jja=(dv(:,2)>=6 & dv(:,2)<=8);
    wbias(i) = nanmean(station{i}.prec(ia(~jja))-station{i}.eiLT(ib(~jja))'.*1e3);
    wrms(i) = sqrt(nanmean((station{i}.prec(ia(~jja))-station{i}.eiLT(ib(~jja))'.*1e3).^2));
    wns(i) = length(cd(~jja))
    wmobs(i) = nanmean(station{i}.prec(ia(~jja)));
    sbias(i) = nanmean(station{i}.prec(ia(jja))-station{i}.eiLT(ib(jja))'.*1e3);
    srms(i) = sqrt(nanmean((station{i}.prec(ia(jja))-station{i}.eiLT(ib(jja))'.*1e3).^2));
    sns(i) = length(cd(jja))
    smobs(i) = nanmean(station{i}.prec(ia(jja)));
%     datestr([station{i}.ei_monthdate(ib(1)),station{i}.ei_monthdate(ib(end))]);
    clear cd ia ib
end


