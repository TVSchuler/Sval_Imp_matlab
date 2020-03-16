%% LT_matlab...
function LT_matlab_2016_09_16


clear all
% close all
clc

% add the other functions to the path...
addpath('./utils/');
addpath('./user/')

disp('Lets get started!!!')

% load user specific variables
[dataset,start,stop,delta_t,datapath,outpath,demname,oro_file,...
    area,region,parallel,nCPU,bgp_thresh,T2file]       = user_specs_LT2014;

% define some constants
g                   = 9.80665;             % gravity acceleration (m s^-2)
rho_w               = 1000;                % density of water in [kg m^-3]
% T0                  = 273.15;              % 0C in K
% R                   = 287;                 % gas constant for dry air [J kg-1 K-1]
% secPerDay           = 86400;               % number of seconds in a day

%% import DEM
disp('read the DEM')
[zzz,X,Y,cellsize,nrow,ncol] = preprocDEM2014(demname);
disp('DEM done')

% the fourier transformed of the terrain
ncx                 = pow2(nextpow2(ncol));
ncy                 = pow2(nextpow2(nrow));
% fft of DEM = hprime; see Smith&Barstad 2004, eq. 49
hprime              = fft2(zzz,ncy,ncx);

%% get the ERA orography
disp('preprocess ERAorography')
[oro,LSoro]                 = preprocORO2014(oro_file,g,area,region,X,Y);
disp('ERAorography done')
% and its fft
oroprime            = fft2(oro,ncy,ncx);

%% make a mask: 200 km around 200m contour of LSoro:
LSmask              = makeLSmask(oro_file,g,area,region);

%% make a time vector
date                = (start:delta_t:stop+1)';
date                = date(date<(stop+1));
% find the months...
dtv                 = datevec(date);
idmon               = find(dtv(:,3)==1&dtv(:,4)==0); %day = 1 and hr =0
% number of months
nmon                = length(idmon);

disp(['going to process ',num2str(nmon),' months'])

% density of air at surface
rho                 = 1.2*exp(-zzz/8500);

% initialize array for monthprec 
monthprec           = zeros(nmon,nrow,ncol);
monthdate           = date(idmon);
monthcorr           = zeros(nmon,nrow,ncol);
%% for each month 
if parallel % parallel computing
    poolobj         = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        parpool(nCPU);
    elseif(poolobj.NumWorkers<nCPU)
        delete(poolobj);
        parpool(nCPU);
    end
    parfor im = 1:nmon
    % for im = 1:nmon
        yyyymm                  = datestr(monthdate(im),'yyyymm');
        disp(['processing month ',yyyymm])
        outfile                 = [outpath,'LT2014_',dataset,'_',region,'_',yyyymm,'.mat'];
%         outfile                 = ['./results_test/LT2015_',dataset,'_',region,'_',yyyymm,'.mat'];
        % check whether this month has been processed before
% % %         if(~exist(outfile,'file'))
            % if not do it
            [precip,time,C]           = LTmonth_fct2015(yyyymm,dataset,datapath,area,g,nrow,ncol,cellsize,zzz,...
                X,Y,bgp_thresh,rho,rho_w,hprime,oroprime,ncx,ncy,delta_t,LSoro,LSmask,T2file);
            mysave(outfile,time,precip);
            % accumulate monthly sum
            monthprec(im,:,:)       = sum(precip,1);
            monthcorr(im,:,:)       = C;
% % %         else
% % %             disp('file exists...omit simulation')
% % %         end
    end % end for each month
% % %     delete(gcp)
else  % sequencial computing
    for im = 1:nmon
        yyyymm                  = datestr(monthdate(im),'yyyymm');
        disp(['processing month ',yyyymm])
        outfile                 = [outpath,'LT2015_',dataset,'_',region,'_',yyyymm,'.mat'];
%         outfile                 = ['./results_test/LT2015_',dataset,'_',region,'_',yyyymm,'.mat'];
        % check whether this month has been processed before
       if(~exist(outfile,'file'))
            % if not do it
            [precip,time,C]           = LTmonth_fct2015(yyyymm,dataset,datapath,area,g,nrow,ncol,cellsize,zzz,...
                X,Y,bgp_thresh,rho,rho_w,hprime,oroprime,ncx,ncy,delta_t,LSoro,LSmask,T2file);
            mysave(outfile,time,precip);
            % accumulate monthly sum
            monthprec(im,:,:)       = sum(precip,1);
            monthcorr(im,:,:)       = C;
        else
            disp('file exists...omit simulation')
        end
    end % end for each month 
end
out2        = ['./results_monthly/LT2015_',datestr(now,30)];
save(out2,'X','Y','monthdate','monthprec','monthcorr')%,'-v7.3')
end % end function

%% embedded functions

% save results for each month
function mysave(outfile,ddd,precip)
    save(outfile,'ddd','precip');
end

% the LTmodel for each month
function [precip,time,C]  = LTmonth_fct2015(yyyymm,dataseries,datapath,area,g,nrow,ncol,cellsize,zzz,...
    X,Y,bgp_thresh,rho,rho_w,hprime,oroprime,ncx,ncy,delta_t,LSoro,LSmask,T2file)
    % seconds in a day, needed to convert units...
    secPerDay           = 86400;
    % read data
    % construct pl_name
    pl_name             = [dataseries,'_pl_',yyyymm];
    % sf_name
    sf_name             = [dataseries,'_sf_',yyyymm];
    % get met data
%     disp('read the met data')
    [T,p,z,RH,U,V,time,lon,lat]     = get_TpzRH_monthly_mat([datapath,pl_name],area,g);
    if(ispc)
        home        = 'M:/';
        datadisk    = 'K:/gggstaff/thomasc/';
    else
        home        = '/uio/kant/geo-natg-u1/thomasc/';
        datadisk    = [home,'datadisk/'];
    end 
    % replace lowes T with T2 from eraInt
    addpath ([home,'matlabtools/snctools2011/mexcdf/mexnc']);
    addpath ([home,'matlabtools/snctools2011/mexcdf/snctools']);
    javaaddpath([home,'matlabtools/snctools2011/netcdfAll-4.2.jar']);
    
   
    ttt             = nc_getall(T2file,Inf);
    ttttime         = datenum(1900,1,1)+double(ttt.time.data)/24;
%     idt             = find(ismember(ttttime,time));
% disp(yyyymm)
    T(:,1,:,:)      = ttt.t2m.data(ismember(ttttime,time),:,:);
    ooo             = repmat(LSoro,1,1,length(time));
    z(:,1,:,:)      = permute(ooo,[3 1 2]);
    clear ttttime ttt idt ooo
    
%     disp('done')
%     disp('calculating the thermodynamics')
    [gamma_m, Gamma_m, Nmsq_m, hw_m, sens_m, q, Tmean, w]     = calculate_thermodyn_4levels_grid_qweighted(T,z,p,RH,g);
%     hw_m = hw_m*2;
%     disp('done')
    % vertical averaging likewise for U and V
%     Umean               = nansum(U.*w,2);
%     Vmean               = nansum(V.*w,2);
%     % index to max windspeed
%     ws                  = sqrt(U.^2+V.^2);
%     [~,iw]              = max(ws(:,2:4,:,:),[],2);
%     Umean               = U(:,iw,:,:);
%     Vmean               = V(:,iw,:,:);
%     Umean               = U(:,2,:,:);
%     Vmean               = V(:,2,:,:);
    Umean               = mean(U(:,2:3,:,:),2);
    Vmean               = mean(V(:,2:3,:,:),2);
    
    mag                 = sqrt(Umean.^2+Vmean.^2);
    % wind direction from 700 hPa...
    dir                 = mod(atan2(U(:,3,:,:),V(:,3,:,:))*180/pi,360);
    
    % re-convert to Umean, Vmean
    Umean               = mag.*sind(dir);
    Vmean               = mag.*cosd(dir);
    
    
    % and for RH
%     RHmean              = nansum(RH.*w,2);
%     RHmean              = RH(:,2,:,:);
    RHmean              = mean(RH(:,2:3,:,:),2);
%     hw_m2               = z(:,4,:,:);
    % clear the averaged variables from workspace
    clear T RH U V z p 

    % calculate tau timescales
    [tauc,tauf]         = calculateTau_lineartransition(hw_m,Tmean);

    % factor to correct for non-saturated conditions, Sinclair(1994)
    RHcorr              = ((RHmean/100-0.8)./0.2).^(0.25); 
    RHcorr(RHmean<=80)  = 0;

    % get the background precip
    [TP,TPdate,TPX,TPY]             = getTP([datapath,sf_name]);
    % check whether both series are synchronized
    if(sum(TPdate-time)>0)
        error('TP and MET not synchronized!')
    end
    
    % number of steps in this month
    nt              = length(time);
    
    % allocate precip matrix
    precip          = zeros(nt,nrow,ncol);

    
    %% for each time step
    for it = 1:nt
        bgpLT           = zeros(size(zzz));
        oprecip         = zeros(size(zzz));
        ttt             = zeros(size(zzz));
        % map TP to the DEM
        % some tests with scatteredInterpolant reveal that it is NOT faster
        % so we stick to griddata (and keep compatibility with older
        % versions of matlab)
        bgp             = griddata(TPX,TPY,squeeze(TP(it,:,:)),X,Y,'cubic'); % in m/6h
        % cut off too light precip
        bgp(bgp<bgp_thresh) = 0;
        
        
% %        horizontal averaging
%         A.Nmsq          = mean(mean(Nmsq_m(it,1,:,:),4),3);
%         A.hw            = mean(mean(hw_m(it,1,:,:),4),3);
%         A.sens          = mean(mean(sens_m(it,1,:,:),4),3);
%         A.tauc          = mean(mean(tauc(it,1,:,:),4),3);
%         A.tauf          = mean(mean(tauf(it,1,:,:),4),3);
%         A.U             = mean(mean(Umean(it,1,:,:),4),3);
%         A.V             = mean(mean(Vmean(it,1,:,:),4),3);
%         A.q             = mean(mean(q(it,1,:,:),4),3);
%         A.rho           = mean(mean(rho)); % rho is only 2D
%         
%         RHfac           = mean(mean(RHcorr(it,1,:,:),4),3);    
%        horizontal averaging, applied only where RHcorr>0

%         ilat = 5:7;ilon=10:25;
%         A.Nmsq          = mean(mean(Nmsq_m(it,1,ilat,ilon),4),3);
%         A.hw            = mean(mean(hw_m(it,1,ilat,ilon),4),3);
%         A.sens          = mean(mean(sens_m(it,1,ilat,ilon),4),3);
%         A.tauc          = mean(mean(tauc(it,1,ilat,ilon),4),3);
%         A.tauf          = mean(mean(tauf(it,1,ilat,ilon),4),3);
%         A.U             = mean(mean(Umean(it,1,ilat,ilon),4),3);
%         A.V             = mean(mean(Vmean(it,1,ilat,ilon),4),3);
%         A.q             = mean(mean(q(it,1,ilat,ilon),4),3);
%         A.rho           = mean(mean(rho(ilat,ilon))); % rho is only 2D
%         RHfac           = mean(mean(RHcorr(it,1,ilat,ilon),4),3);          
%         
        RHfac           = squeeze(RHcorr(it,1,:,:));   
%         RHmask          = double(LSoro>5);
% make a mask for averaging: 200km around the 200m-contour + where
% TP>bgp_thresh
%         mask = double(squeeze(TP(it,:,:))>bgp_thresh).*LSmask;
% % %         mask = zeros(size(LSmask));
% % %         mask = double(squeeze(TP(it,:,:))>0).*LSmask;
        mask = LSmask;
%         mask = LSmask;
        mask(mask==0) = NaN;
%          RHmask(LSoro>200) = NaN;
%          RHmask(LSoro<5) = NaN;
%         RHmask(RHfac==0) = NaN;
%         bgptmp=squeeze(TP(it,:,:));
%         RHmask(bgptmp==0) = NaN;
%         RHmask          = double(RHfac>0);
%   RHmask1 = NaN(size(RHfac));       
% 
% m_u   = median(median(Umean(it,1,:,:))); m_v = median(median(Vmean(it,1,:,:)));        
% m_dir = mod(atan2(-m_u,-m_v)*180/pi,360); % dominating wind dir (in meteo sense, i.e. where it is coming from)
%     % output from atan2d is in degrees [-180, 180], convert to [0 360]
%     % and turn it by 22.5 deg such that limits of classes are at 
%     % 22.5, 67.5...and not at 45, 90 etc
%     % atan2d was introduced in R2012b, to keep compatibility with older
%     % versions, we use atan2 and do the rad2deg conversion
% 
%     % classify into 8 classes:
%     m_dir          = floor((m_dir+22.5)/45);
%     m_dir          = mod(m_dir,8)+1; % 337.5-22.5 = 1, N....classes from 1-8
%     
%     switch m_dir
%         case 1 % N
%             RHmask1(1:5,10:30) = 1;
%         case 2 % NE
%             RHmask1(1:5,20:end) = 1;
%         case 3 % E
%             RHmask1(3:8,20:end) = 1;
%         case 4 % SE
%             RHmask1(5:end,20:end) = 1;
%         case 5 % S
%             RHmask1(5:end,10:30) = 1;
%         case 6 % SW
%             RHmask1(5:end,1:20) = 1;
%         case 7 % W
%             RHmask1(3:8,1:20) = 1;
%         case 8 % NW
%             RHmask1(1:5,1:20) = 1;
%         otherwise
%             error('tull med vind')
%     end
%     
%         
% RHmask = RHmask.*RHmask1;

        RHfac           = RHfac.*mask;
        RHfac           = nanmean(RHfac(:));

        A.Nmsq          = squeeze(Nmsq_m(it,1,:,:)).*mask;
        A.hw            = squeeze(hw_m(it,1,:,:)).*mask;
        A.sens          = squeeze(sens_m(it,1,:,:)).*mask;
        A.tauc          = squeeze(tauc(it,1,:,:)).*mask;
        A.tauf          = squeeze(tauf(it,1,:,:)).*mask;
        A.U             = squeeze(Umean(it,1,:,:));%.*mask;
        A.V             = squeeze(Vmean(it,1,:,:));%.*mask;
        A.q             = squeeze(q(it,1,:,:)).*mask;
        A.rho           = rho; % rho is only 2D
       
        A.Nmsq          = nanmean(A.Nmsq(:));
        A.hw            = nanmean(A.hw(:));
        A.sens          = nanmean(A.sens(:));
        A.tauc          = nanmean(A.tauc(:));
        A.tauf          = nanmean(A.tauf(:));
        A.U             = nanmean(A.U(:));
        A.V             = nanmean(A.V(:));
        A.q             = nanmean(A.q(:));
        A.rho           = nanmean(A.rho(:));

        % bgp correction
        if(RHfac>0)
            % apply to ERAorography
            p1          = apply_LTmodel2014(A,oroprime,ncx,ncy,ncol,nrow,cellsize);
            bgpLT       = p1.*(delta_t*secPerDay).*RHfac/rho_w; % convert kg/s/m2 to m/6h
%             % correct for drying of the air (Smith&Evans, 2009)
            WV          = vaporflux(A.rho,A.q,A.U,A.V,A.hw);
            DR1          = drying_ratio(p1,WV,A.U,A.V,cellsize);
            p1pos        = p1>0;
            ttt(p1pos)   = bgpLT(p1pos).*(1-exp(-DR1(p1pos)));
            bgpLT(p1pos) = bgpLT(p1pos).*exp(-DR1(p1pos));

            % apply LT to 1km DEM
            p2          = apply_LTmodel2014(A,hprime,ncx,ncy,ncol,nrow,cellsize);
            % apply the mask
            oprecip     = p2.*(delta_t*secPerDay).*RHfac/rho_w;
            % correct for drying of the air (Smith&Evans, 2009)
%             WV          = vaporflux(A.rho,A.q,A.U,A.V,A.hw);
            DR          = drying_ratio(p2,WV,A.U,A.V,cellsize);
            % correct for drying
            ppos              = p2>0;
            ttt(ppos)         = (1-exp(-DR(ppos))).*oprecip(ppos);
            oprecip(ppos)     = oprecip(ppos).*exp(-DR(ppos));
        end

        % correct bgp and add oprecip
        corrbgp             = bgp;% - bgpLT;
        corrbgp             = max(corrbgp,0);           % do not allow negative values
        tmpprecip           = corrbgp +oprecip;
        tmpprecip           = max(tmpprecip,0);         % do not allow negative values
        precip(it,:,:)      = tmpprecip;
        if it==1
            C = ttt;
        else
            C = C+ttt;%mean(cat(3,C,ttt),3);
        end
    end % end for each time step
   

end