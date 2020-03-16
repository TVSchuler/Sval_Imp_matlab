function downscaleERA_matlab_TFiddes_RHwindInterp_addcorrection
% patch to correct for unrealistically high temperature when T>0C:
% we use 2D-interpolation of ERA T2 instead.

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
T0          = 273.15;
%% define time period
dataset  = 'ei';
switch dataset
    case 'ei'

        start       = datenum(1979,1,1); % ei data from 1.1.1979
        stop        = datenum(2018,12,31); % until 31.3.2013
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
datapathT2      = [datadisk,'eraSvalbard/ecmwf-api-client-python/'];

% path for output files
outpath     = [datadisk2,'downscaledERA_PL_SRF/',dataset,'/'];

% define area of interest: [min(lon) max(lon) max(lat) min(lat)]
area        = [6 36 75 82.5]; % entire Svalbard according to ERAdata

%% define filenames (USER SPECIFIC)
demname             = [datadisk,'austfonna/dem_sval1km/sval1km_wgs84.txt'];

%----------------------------------------------------------------------------
% ERAint orography
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
month_date          = date(idmon);

%% open the T2 file (hard-wired file: USER adjust to your case)
    ttt             = nc_getall([datapath,dataset,'_Svalbard_T2.nc'],Inf);
 % open the land-sea mask file 
    mmm             = nc_getall([datapathT2,'data/',dataset,'_landsea_mask.nc'],Inf);
   lsm = mmm.lsm.data;
% make the LScoordinates   
    [LON,LAT]        = meshgrid(lon,lat);

    % convert to UTM33
    [LSX,LSY,~]      = deg2utm33(LAT(:),LON(:));
    LSX              = reshape(LSX,length(lat),length(lon));
    LSY              = reshape(LSY,length(lat),length(lon));

%% for each month % parfor?
% % s           = matlabpool('size');
% % if(s==0)
% %   matlabpool('open',12)
% % end
% % parfor im = 1:nmon
for im = 1:nmon
    yyyymm                  = datestr(month_date(im),'yyyymm');
    disp(['processing month ',yyyymm])
    outfile                 = [outpath,'ERA_Tpl_RHsrf_UVsrf_',dataset,'_Svalbard_6h_',yyyymm,'.mat'];
    ERA_correctTfiddes(outfile,ttt,X,Y,LSX,LSY,T0,lsm);
end
    
end
function ERA_correctTfiddes(outfile,ttt,X,Y,LSX,LSY,T0,lsm)
        srf_date = double(ttt.time.data)/24+datenum(1900,1,1);
        load(outfile,'ddd','temp'); % ddd, temp, (rh, u, v)
        tflag                   = 0; % flag to indicate whether temp was updated or not
        % intersect the dates
        [dd2,ia,ib]         = intersect(ddd,srf_date);
        if(isempty(dd2))
            error('error matching dates')
        end    
        for it=1:length(dd2) % current step of ddd is ia(it) and of srf_date: ib(it)
            % check whether current temp>0?
            idx = find(temp(ia(it),:,:)>=T0);
            if(~isempty(idx))
                Ttmp             = squeeze(temp(ia(it),:,:));
                F                = scatteredInterpolant(LSX(:),LSY(:),reshape(ttt.t2m.data(ib(it),:,:),[],1));
                T2int            = F({X(1,:),Y(:,1)})';
                % fill in
                Ttmp(idx)        = T2int(idx);
                temp(ia(it),:,:) = Ttmp;
                tflag            = 1;
            end   
        end
        if(tflag) % only save if we have made changes...
            disp(['updating ',outfile])
            save(outfile,'ddd','temp','-append') % append --> overwrite the changes 
        end
end