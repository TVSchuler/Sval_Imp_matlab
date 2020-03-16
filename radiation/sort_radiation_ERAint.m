% sort_radiation_ERA
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
% addpath('./utils/');
% enable SNCTOOLS & mexnc
addpath ([home,'matlabtools/snctools2011/mexcdf/mexnc']);
addpath ([home,'matlabtools/snctools2011/mexcdf/snctools']);
javaaddpath([home,'matlabtools/snctools2011/netcdfAll-4.2.jar']);

% for ERAint
% path for output files
dataset     = 'ei';
outpath     = [datadisk2,'downscaledERA_RAD/',dataset,'/'];

%% hard-wired path...requires USER adjustment!

% %% load data: SWtoa
tmp   = nc_getall([datadisk,'eraSvalbard/eraInt/ei_Svalbard_TOA.nc'],Inf);
date  = datenum(1900,1,1)+(double(tmp.time.data)/24);

% pick out every second value to get 6h data
date_6h  = date(8:2:end);

SWtoa_6h = tmp.tisr.data(8:2:end,:,:);
SWtoa_6h(1,:,:) = SWtoa_6h(1,:,:)-tmp.tisr.data(6,:,:);
% do differencing to convert cumulative to actual values....
SWtoa_6h(3:2:end,:,:) = SWtoa_6h(3:2:end,:,:)-SWtoa_6h(2:2:end-1,:,:);
% convert from J/m2 to W/m2
SWtoa_6h      = SWtoa_6h/(6*3600);

%% hard-wired path...requires USER adjustment!
% load data: SWsfc
tmp   = nc_getall([datadisk,'eraSvalbard/eraInt/ei_Svalbard_SWin2018_3.nc'],Inf);
dat2  = datenum(1900,1,1)+(double(tmp.time.data)/24);

% pick out every second value to get 6h data
dat2_6h  = dat2(8:2:end);
date_6h = dat2_6h;

SWsfc_6h = tmp.ssrd.data(8:2:end,:,:);
SWsfc_6h(1,:,:) = SWsfc_6h(1,:,:)-tmp.ssrd.data(6,:,:);
% do differencing to convert cumulative to actual values....
SWsfc_6h(3:2:end,:,:) = SWsfc_6h(3:2:end,:,:)-SWsfc_6h(2:2:end,:,:);
% convert from J/m2 to W/m2
SWsfc_6h      = SWsfc_6h/(6*3600);

%% hard-wired path...requires USER adjustment!
% load data: LWsfc
tmp   = nc_getall([datadisk,'eraSvalbard/eraInt/ei_Svalbard_LWin2018_3.nc'],Inf);
dat3  = datenum(1900,1,1)+(double(tmp.time.data)/24);

% pick out every second value to get 6h data
dat3_6h  = dat3(8:2:end);

LWsfc_6h = tmp.strd.data(8:2:end,:,:);
LWsfc_6h(1,:,:) = LWsfc_6h(1,:,:)-tmp.strd.data(6,:,:);

% do differencing to convert cumulative to actual values....
LWsfc_6h(3:2:end,:,:) = LWsfc_6h(3:2:end,:,:)-LWsfc_6h(2:2:end,:,:);
% convert from J/m2 to W/m2
LWsfc_6h      = LWsfc_6h/(6*3600);

% write all into monthly files
% find the months...
dtv                 = datevec(date_6h);
dvm                 = dtv(:,1)*100+dtv(:,2);
udvm                = unique(dvm);

nmon                = length(udvm);

for im = 1:nmon-1  % we do not process the only value for Jan 2014...
% for im = 1:nmon
        outfile     = [outpath,'ERA_RAD_',dataset,'_Svalbard_6h_',num2str(udvm(im)),'.mat'];
        idm         = find(dvm==udvm(im));
%         filename    = [datapath,dataset,'_pl_',num2str(udvm(im))];
        era         = struct('longitude',tmp.longitude.data,'latitude',tmp.latitude.data,...
            'time',date_6h(idm),'SWtoa_6h',SWtoa_6h(idm,:,:),'SWsfc_6h',SWsfc_6h(idm,:,:),'LWsfc_6h',LWsfc_6h(idm,:,:));
        %             'time',date_6h(idm),'SWsfc_6h',SWsfc_6h(idm,:,:),'LWsfc_6h',LWsfc_6h(idm,:,:));

        save(outfile,'era')
        disp(['saved to ',outfile])

end