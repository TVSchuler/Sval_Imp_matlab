function [dataset,start,stop,delta_t,datapath,outpath,...
    demname,oro_file,area,region,parallel,nCPU,bgp_thresh,T2file] ...
    = user_specs_LT2014
% user specifications for LT_matlab_2014
if(ispc)
    home        = 'M:/';
%     datadisk    = 'K:/gggstaff/thomasc/';
    datadisk    = 'L:/';
%     datadisk2    = 'K:/gggruppe-prosj/glasio-data/';
    datadisk2   = 'V:/data/austfonna/processed/';
else
    home        = '/uio/kant/geo-natg-u1/thomasc/';
    datadisk    = [home,'datadisk/'];
    datadisk2   = [home,'glacier_data/'];
end

%% define dataset, time period and path to ERA data
%  'ei' = ERAinterim, 'e4' = ERA40
dataset  = 'ei';

switch dataset
    case 'ei'
%         start       = datenum(1979,1,1); % ei data from 1.1.1979
% %         stop        = datenum(2013,3,31); % until 31.3.2013
%          start       = datenum(2015,7,1); % ei data from 1.1.1979
%         stop        = datenum(2015,6,30); % until 31.3.2013
%         stop        = datenum(2016,6,30); % until 31.3.2013
        start       = datenum(2018,9,1); % ei data from 1.1.1979
        stop        = datenum(2018,12,31); % until 31.3.2013
%         start       = datenum(2003,9,1);
%         stop        = datenum(2011,4,30);
%        stop        = datenum(1979,3,31); % test
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

% path for output files
% outpath     = [datadisk2,'svalbard/austfonna/precipLT2015/',dataset,'/'];
outpath     = [datadisk2,'precipLT2015/',dataset,'/'];


%% define filenames DEM
% (entire path, in case the DEM is not at the same location as ERAdata)
% the DEM should be in ARC-ascii raster format
% coordinates should be UTM, otherwise some modifications in the code are
% required (several instances, each time when converting lat/lon to the
% grid-coordinates...
demname             = [datadisk,'austfonna/dem_sval1km/sval1km_wgs84.txt'];
%----------------------------------------------------------------------------
% ERAint orography
oro_file            = 'ERAintOrography2014'; % located in the ./utils-folder
% define area of interest: [min(lon) max(lon) max(lat) min(lat)]
area        = [6 36 75 82.5]; % entire Svalbard 
region      = 'svalbard';     % ID to be used in output file name and for spec functions (if applicable)

%% multicore option: if matlab parallel-computing toolbox available and computer has several CPUs
parallel    = false;  % specify: true OR false
nCPU        = 6;     % number of CPUs to be used

%% threshold for TP correction; 
% all TP below that value is discarded, was used previously (ERA40 in Iceland had 
% a problem with overestimating low intensity precip)
% % bgp_thresh          = 0.0017; % m/delta_t
bgp_thresh          = 0.00015; % m/delta_t

T2file     = [datadisk,'eraSvalbard/eraInt/ei_Svalbard_T2_2018_3.nc'];
end % end of function, list of variables can be extended but must be specified in the output
