function downscaleERA_matlab_TFiddes_RHwindInterp
% downscale ERA fields to high-res DEM following TopoSCALE by Fiddes&Gruber
% (2014): for temperature
% modification: in case T>0C: interpolation of ERA field to high-res DEM
% RH and windspeed is also interpolated without altitude adjustment

clear all
close all
clc

% USER SPECIFIC: adjust to your case
if(ispc)
    home        = 'M:/';
    datadisk    = 'L:/';
    datadisk2    = 'V:/data/example/processed/';
else
    home        = '/SERVER/DISK/USER/HOME/';
    datadisk    = [home,'datadisk/'];
    datadisk2   = [home,'glacier_data/'];
end
% enable asciigridread etc...
addpath(strcat(home,'matlabtools/'));          %% 
addpath('./utils/');
% enable SNCTOOLS & mexnc
addpath ([home,'matlabtools/snctools2011/mexcdf/mexnc']);
addpath ([home,'matlabtools/snctools2011/mexcdf/snctools']);
javaaddpath([home,'matlabtools/snctools2011/netcdfAll-4.2.jar']);

disp('Lets get started!!!')

% define some constants
g           = 9.80665;             % gravity acceleration (m s^-2)

%% define time period
dataset  = 'ei';
switch dataset
    case 'ei'
        start       = datenum(1979,1,1); % ei data from 1.1.1979
        stop        = datenum(2018,12,31); % until 31.12.2018
        delta_t     = 6/24; % 6h temp resolution of ERA
        datapath    = [datadisk,'eraSvalbard/eraInt/'];
    case 'e4'
        start       = datenum(1957,9,1); % e40 data from 1.9.1957
        stop        = datenum(2002,8,31); % until 31.8.2002
        delta_t     = 6/24; % 6h temp resolution of ERA
        datapath    = [datadisk,'eraSvalbard/era40/'];
    otherwise
        error('dataseries not defined!')
end

% path to srf-ncfiles (USER SPECIFIC)
datapathT2      = [datadisk,'eraSvalbard/eraInt/'];

% path for output files
outpath     = [datadisk2,'downscaledERA_PL_SRF/',dataset,'/'];

% define area of interest: [min(lon) max(lon) max(lat) min(lat)]
area        = [6 36 75 82.5]; % entire Svalbard according to ERAdata

%% define filenames (USER SPECIFIC)
demname             = [datadisk,'austfonna/dem_sval1km/sval1km_wgs84.txt'];
% demname             = [datadisk,'austfonna/dem_sval1km/na_1km_wgs84.txt'];

%----------------------------------------------------------------------------
% ERAint orography (USER SPECIFIC)
datapath1           = [datadisk,'austfonna/eraInt/download/data/'];
oro_file            = [datapath1,'ei_orography'];

%% import DEM
disp('read the DEM')
[zzz,X,Y,~,nrow,ncol,~,~] = preprocDEM(demname);
disp('DEM done')

%% get the lat-lon-level from the era-files:
% pick one for a date that exists in both...
tmpdate                  = '200001';
tmpname                  = [datapath,dataset,'_pl_',tmpdate];

[lat,lon,~]              = getERAtmp(tmpname,area);


%% get the ERA orography
disp('preprocess ERAorography')
oro                 = preprocORO_latlon(oro_file,g,area,lat,lon);
disp('ERAorography done')

%% make a time vector
date                = (start:delta_t:stop+1)';
date                = date(date<(stop+1));
% find the months...
dtv                 = datevec(date);
idmon               = find(dtv(:,3)==1&dtv(:,4)==0); %day = 1 and hr =0
% number of months
nmon                = length(idmon);

% initialize array for monthprec % do we still need this?
month_temp          = zeros(nmon,nrow,ncol);
month_rh            = zeros(nmon,nrow,ncol);
month_u             = zeros(nmon,nrow,ncol);
month_v             = zeros(nmon,nrow,ncol);
month_date          = date(idmon);

 % open the T2 file 
    ttt             = nc_getall([datapathT2,dataset,'_Svalbard_T2_2018_3.nc'],Inf);
    rrr             = nc_getall([datapathT2,dataset,'_Svalbard_T2dew2018_3.nc'],Inf);
    www             = nc_getall([datapathT2,dataset,'_Svalbard_wind_10m_2018_3.nc'],Inf);
    
    % convert T2d to vapour pressure: Bolton 1980
    e_s0   = 6.112;                           % e_s at T=0 C [hPa];
    T_0    = 273.15; % K2degC
    pvap =    e_s0*exp(17.67*(rrr.d2m.data-T_0)./(rrr.d2m.data-29.65));
    [ddd,ia,ib] = intersect(ttt.time.data,rrr.time.data);
    
    psat =    e_s0*exp(17.67*(ttt.t2m.data(ia,:,:)-T_0)./(ttt.t2m.data(ia,:,:)-29.65));
    
    relHum   = pvap./psat;
    u10  = www.u10.data;
    v10  = www.v10.data;
    
    srf_date = double(ttt.time.data(ia))/24+datenum(1900,1,1);
    srf_dat2 = double(www.time.data)/24+datenum(1900,1,1);
    



%% for each month % parfor?
nCPU = 12;
    poolobj         = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
%         parpool(nCPU);
    elseif(poolobj.NumWorkers<nCPU)
        delete(poolobj);
        parpool(nCPU);
    end

% parfor im = 1:nmon
 for im = 1:nmon
    yyyymm                  = datestr(month_date(im),'yyyymm');
    disp(['processing month ',yyyymm])
    outfile                 = [outpath,'ERA_Tpl_RHsrf_UVsrf_',dataset,'_Svalbard_6h_',yyyymm,'.mat'];
%     outfile                 = ['LTcluster_',dataset,'_Svalbard_6h_',yyyymm,'.mat'];
    % check whether this month has been processed before
    if(~exist(outfile,'file'))
        % if not do it
        % extract the actual month from the SRF data
        if im<nmon
            idx = find(srf_date>=month_date(im)&srf_date<month_date(im+1));
            id2 = find(srf_dat2>=month_date(im)&srf_dat2<month_date(im+1));
        elseif im==nmon
            idx = find(srf_date>=month_date(im));
            id2 = find(srf_dat2>=month_date(im));
        else
            error('something went wrong with matching SRFdate and monthdate')
        end
        
        RHsrf                   = relHum(idx,:,:);
        Usrf                    = u10(id2,:,:);
        Vsrf                    = v10(id2,:,:);
        
        [temp,rh,u,v,time]           = ERA_PL_SRF_month_fct(yyyymm,dataset,datapath,area,g,nrow,ncol,zzz,X,Y,RHsrf,Usrf,Vsrf);
        % check whether times coincide...
%         ddd = intersect(time,srf_date(idx));
%         if length(ddd)~=length(idx)
%             error('problem matching dates...')
%         end
        mysave(outfile,time,temp,rh,u,v);
        % accumulate monthly sum
        month_temp(im,:,:)       = mean(temp,1);
        month_rh(im,:,:)         = mean(rh,1);
        month_u(im,:,:)          = mean(u,1);
        month_v(im,:,:)          = mean(v,1);
    else
        disp('file exists...omit simulation')
    end
end % end for each month
% matlabpool('close')
out2        = ['ERA_PLds_Fiddes_',datestr(now,30)];
save(out2,'X','Y','month_date','month_temp','month_rh','month_u','month_v','-v7.3')
end % end function
function mysave(outfile,ddd,temp,rh,u,v)
    save(outfile,'ddd','temp','rh','u','v');
end
function [temp,rh,u,v,time]  = ERA_PL_SRF_month_fct(yyyymm,dataset,datapath,area,g,nrow,ncol,zzz,X,Y,RHsrf,Usrf,Vsrf);

    % read data
    % construct pl_name
    pl_name         = [dataset,'_pl_',yyyymm];
 
    % get met data
    [T,~,z,~,~,~,time,lon,lat]   = get_TpzRH_monthly_mat([datapath,pl_name],area,g);
    % number of vertical levels
    sT               = size(T);
    nlev             = sT(2);
    nt               = sT(1);
%     nlon             = sT(4);
%     nlat             = sT(3);
    
    % make the LScoordinates   
    [LON,LAT]        = meshgrid(lon,lat);

    % convert to UTM33
    [LSX,LSY,~]      = deg2utm33(LAT(:),LON(:));
    LSX              = reshape(LSX,length(lat),length(lon));
    LSY              = reshape(LSY,length(lat),length(lon));

    % prepare input for griddatan: matrices must be nlev+1*nlat*nlon
    % constant vars:
    LSX2             = repmat(LSX,[1 1 nlev+1]); %LSX2 = permute(LSX2,[3 1 2]); % lev lat lon
    LSY2             = repmat(LSY,[1 1 nlev+1]); %LSY2 = permute(LSY2,[3 1 2]);
    
    % allocate memory for result Ts
    temp             = zeros(nt,nrow,ncol);
    rh               = zeros(nt,nrow,ncol);
    u                = zeros(nt,nrow,ncol);
    v                = zeros(nt,nrow,ncol);
    
    % time loop
    
    for it = 1:nt
        % inside the time loop
        % consider actual time step
        z2               = squeeze(z(it,:,:,:)); 
        T2               = squeeze(T(it,:,:,:)); 
%         RH2              = squeeze(RH(it,:,:,:)); 
%         U2               = squeeze(U(it,:,:,:)); 
%         V2               = squeeze(V(it,:,:,:)); 
        % permute 
        z2               = permute(z2,[2 3 1]);
        T2               = permute(T2,[2 3 1]);
%         RH2              = permute(RH2,[2 3 1]);
%         U2               = permute(U2,[2 3 1]);
%         V2               = permute(V2,[2 3 1]);
        % add an extra layer at bottom to avoid extrapolation problem
        z2               = cat(3,zeros(size(LSX))+min(zzz(:)),z2); % add min of DEM
        T2               = cat(3,T2(:,:,1),T2); % repeat lowermost layer
%         RH2              = cat(3,RH2(:,:,1),RH2); % repeat lowermost layer
%         U2               = cat(3,U2(:,:,1),U2); % repeat lowermost layer
%         V2               = cat(3,V2(:,:,1),V2); % repeat lowermost layer

%         T_tmp            = griddatan([LSX2(:),LSY2(:),z2(:)],T2(:),[X(:),Y(:),zzz(:)]);
        F                = scatteredInterpolant([LSX2(:),LSY2(:),z2(:)],T2(:));
        F2               = scatteredInterpolant(LSX(:),LSY(:),reshape(RHsrf(it,:,:),[],1));
        F3               = scatteredInterpolant(LSX(:),LSY(:),reshape(Usrf(it,:,:),[],1));
        F4               = scatteredInterpolant(LSX(:),LSY(:),reshape(Vsrf(it,:,:),[],1));

%         F2               = scatteredInterpolant([LSX2(:),LSY2(:),z2(:)],RH2(:));
%         F3               = scatteredInterpolant([LSX2(:),LSY2(:),z2(:)],U2(:));
%         F4               = scatteredInterpolant([LSX2(:),LSY2(:),z2(:)],V2(:));
        % evaluate it
        T_tmp            = F([X(:),Y(:),zzz(:)]);
        RH_tmp           = F2(X(:),Y(:));
        U_tmp            = F3(X(:),Y(:));
        V_tmp            = F4(X(:),Y(:));
        % reshape it to the grid and put it into matrix
        temp(it,:,:)     = reshape(T_tmp,size(zzz));
        rh(it,:,:)       = reshape(RH_tmp,size(zzz));
        u(it,:,:)        = reshape(U_tmp,size(zzz));
        v(it,:,:)        = reshape(V_tmp,size(zzz));
        
    end
    

    

end
