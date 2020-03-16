% sort_T2_T2dew_ERA
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
outpath     = [datadisk2,'ERA_Svalbard/',dataset,'/'];

%% hard-wired path...requires USER adjustment!
% load data: ERA-T2
tmp   = nc_getall([datadisk,'eraSvalbard/eraInt/ei_Svalbard_T2.nc'],Inf);
date  = datenum(1900,1,1)+(double(tmp.time.data)/24);

T2  = tmp.t2m.data;

%% hard-wired path...requires USER adjustment!
% load data: ERA-T2dewpoint
tmp   = nc_getall([datadisk,'eraSvalbard/eraInt/ei_Svalbard_T2dew.nc'],Inf);
dat2  = datenum(1900,1,1)+(double(tmp.time.data)/24);

T2dew = tmp.d2m.data;

% convert T2d to vapour pressure: Bolton 1980
e_s0   = 6.112;                           % e_s at T=0 C [hPa];
T_0    = 273.15; % K2degC
pvap =    e_s0*exp(17.67*(T2dew-T_0)./(T2dew-29.65));

% write all into monthly files
% find the months...
dtv                 = datevec(date);
dvm                 = dtv(:,1)*100+dtv(:,2);
udvm                = unique(dvm);

nmon                = length(udvm);

% for im = 1:nmon-1  % we do not process the only value for Jan 2014...
for im = 1:nmon
        disp(['processing month: ',num2str(udvm(im))])
        outfile     = [outpath,'ERA_T2_pvap_',dataset,'_Svalbard_6h_',num2str(udvm(im)),'.mat'];
        if(~exist(outfile,'file'))
            idm         = find(dvm==udvm(im));
    %         filename    = [datapath,dataset,'_pl_',num2str(udvm(im))];
            era         = struct('longitude',tmp.longitude.data,'latitude',tmp.latitude.data,...
                'time',date(idm),'T2',T2(idm,:,:),'pvap',pvap(idm,:,:));
            save(outfile,'era')
            disp(['saved to ',outfile])
        else
            disp('file exists...omit simulation')
        end

end