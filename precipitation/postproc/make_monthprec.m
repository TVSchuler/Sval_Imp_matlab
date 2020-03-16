% make monthprec by summing precip fields for each month
clear all
close all
clc

if(ispc)
    home        = 'M:/';
    datadisk    = 'K:/gggstaff/thomasc/';
    datadisk2    = 'V:/data/austfonna/processed/precipLT2015/';
else
    home        = '/uio/kant/geo-natg-u1/thomasc/';
    datadisk    = [home,'datadisk/'];
    datadisk2   = [home,'glacier_data/'];
end
% enable asciigridread etc...
addpath(strcat(home,'matlabtools/'));          %% 

for id = 1

%% define time period
if(id==1)
    dataset  = 'ei';
elseif(id==2)
    dataset = 'e4';
else
    error('undefined dataset')
end
disp(['processing dataset ',dataset])
switch dataset
    case 'ei'
        start       = datenum(1979,1,1); % ei data from 1.1.1979
%         stop        = datenum(2015,6,30); % until 31.3.2013
%         start       = datenum(2015,7,1); % ei data from 1.1.1979
        stop        = datenum(2018,8,31); % until 31.3.2013
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
outpath     = [datadisk2,dataset,'/'];
%% make a time vector
date                = (start:delta_t:stop+1)';
date                = date(date<(stop+1));
% find the months...
dtv                 = datevec(date);
idmon               = find(dtv(:,3)==1&dtv(:,4)==0); %day = 1 and hr =0
% number of months
nmon                = length(idmon);
% preallocate...
monthprec           = zeros(nmon,548,448); % size of the DEM * nmon
monthdate           = date(idmon);

% for each month
for im = 1:nmon
    yyyymm                  = datestr(monthdate(im),'yyyymm');
    disp(['processing month ',yyyymm])
    outfile                 = [outpath,'LT2014_',dataset,'_svalbard_',yyyymm,'.mat'];
    % check whether this month has been processed before
    if(exist(outfile,'file'))
        load(outfile); % time, precip
        % accumulate monthly sum
        monthprec(im,:,:)       = sum(precip,1);
    else
        disp('file does not exist...omitting month')
    end
end

save(['./results_monthly/LT2015_',dataset,'_Svalbard_monthly_19792018'],'monthdate','monthprec')

end