% read ERA-fields and correct SW according to Fiddes & Gruber 2014
%
% 1) partion SW_grid into SWdir and SWdif
% 2) elevation adjustment of SWdir, dz between ERA-grid and Sval1km
% 3) Topographic correction of SWdir and SWdif
%

clear all
close all

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
%% paths (USER-SPECIFIC: adjust to your case)
addpath('./m_map');
ERARAD_path = [datadisk2,'downscaledERA_RAD/ei/ERA_RAD_ei_Svalbard_6h_'];
ERAsub_path = [datadisk2,'downscaledERA_PL_SRF/ei/ERA_Tpl_RHsrf_UVsrf_ei_Svalbard_6h_'];
ERAgrid_path = [datadisk2,'ERA_Svalbard/ei/ERA_T2_pvap_ei_Svalbard_6h_'];
sun_path = [datadisk2,'downscaledSUN/SolarGeometry_6h_MOMENTAN/SolarGeometry_jd'];
out_path = [datadisk2,'downscaledSUN/downscaled_ERAi_SW_LW_6h_momentan/SWLW_'];

%% settings
% t_beg = '201506';%'201311';%'197902';
% t_beg = '197901'; % january79 is missing the first val...
% t_end = '201606';%'201303';
t_beg = '201811';
t_end = '201812';
P0 = 101325;%  	sea level standard atmospheric pressure 	101325 Pa
%L 	temperature lapse rate, = g/cp for dry air 	0.0065 K/m
%cp 	constant pressure specific heat 	~ 1007 J/(kg•K)
T0 = 288.15;%	sea level standard temperature 	288.15 K
g = 9.80665;% 	Earth-surface gravitational acceleration 	9.80665 m/s2
M = 0.0289644;%	molar mass of dry air 	0.0289644 kg/mol
R = 8.31447; %	universal gas constant 	8.31447 J/(mol•K)
aa = 0.43; % emissivity paramterization coef  Gubler et al 2012
bb = 5.7; % emissivity paramterization coef, Gubler et 2012

sigma = 5.6704e-8; % stafan-boltzman
Rv = 461.5; % J K /kg  Gas konstant for water vapor
Lv = 2.5e6; %J/kg Latent heat of evaporation
To = 273.15; %K, ref Temp
es0=611; % Pa, ref saturation vapor pressure

infoRAD= {'Downscaled ERAi SW and LW to Svalbard 1000m following Fiddes & Gruber 2014';'created by "ERA_SW_to_Sval1km.m" using solar geomtry from SolarGeomtry.m';'6h avg of the last 6 hours';'No data on edges for SW';'Zenith > 89.5 SW==0'};

%% load
% the high resolution DEM
load([datadisk2,'downscaledSUN/data/Sval_DEM1000m.mat'])
% the ERA-DEM 
load([datadisk2,'downscaledSUN/data/ERAintOrography.mat']);

%% pre-proc GRIDs and TIME
% compute difference between ERA and Sval1km topo
demERA = interpERA2SVAL(LSlon,LSlat,oro2014,x,y);
dz=dem-demERA;
ncols = size(dem,2);
nrows = size(dem,1);

[eraLON,eraLAT] = meshgrid(double(LSlon),double(LSlat));


% pressure
Pera = P0 .* exp(-(g*M*demERA)/(R*T0));
Psub = P0 .* exp(-(g*M*dem)/(R*T0));
dP = Psub-Pera;



% clear all redundant info
clear Pera Psub dem demERA ni nj nk aspect FGC slope
close all

% set time vector
time = datenum(t_beg,'yyyymm'):1:datenum(t_end,'yyyymm');
[tmp1,tmp2,~,~,~,~]=datevec(time);
time=unique(datenum(tmp1,tmp2,ones(size(tmp1))));

%break

% % %% loop over months
% poolobj         = gcp('nocreate'); % If no pool, do not create new one.
%     if isempty(poolobj)
%         poolobj=parpool(6);
%     elseif(poolobj.NumWorkers<6)
%         delete(poolobj);
%         poolobj=parpool(6);
%     end
% % pool=parpool(12);
% % matlabpool open
% parfor it=1:length(time)
for it=1:length(time)
%i=1;
% for it=1:length(time)
    t=time(it);
%     disp(datestr(t))
    
    if(~exist([out_path,datestr(t,'yyyymm'),'.mat'],'file'))% && i~=173) % only if not already done
        disp(datestr(t))
        % load monthly downscaled ERA-data (sub)
        ERAsub=load([ERAsub_path,datestr(t,'yyyymm'),'.mat']);
        % load orginal ERA-data (grid)
        ERAgridrad=load([ERARAD_path,datestr(t,'yyyymm'),'.mat']);
        ERAgridTvP=load([ERAgrid_path,datestr(t,'yyyymm'),'.mat']);
        
        % time syncron test
        if(sum(abs(diff([length(ERAsub.ddd),length(ERAgridrad.era.time),length(ERAgridTvP.era.time)])))~= 0)
            %-size(ERAsub.ddd)+size(ERAgridrad.era.time)-size(ERAgridTvP.era.time)) ~= 0)
            disp('-------Time not syncron------'),
            disp(datestr(t))
            if(t==datenum(1979,1,1)||t==datenum(2013,4,1))
                % skip first data
                ERAsub.ddd = ERAsub.ddd(2:end);
                ERAsub.temp = ERAsub.temp(2:end,:,:);
                ERAsub.rh = ERAsub.rh(2:end,:,:);
                ERAgridTvP.era.time = ERAgridTvP.era.time(2:end);
                ERAgridTvP.era.T2 = ERAgridTvP.era.T2(2:end,:,:);
                ERAgridTvP.era.pvap = ERAgridTvP.era.pvap(2:end,:,:);
            end
            if(t==datenum(2016,7,1)||t==datenum(2017,1,1)||t==datenum(2018,1,1)||t==datenum(2018,9,1))
                % skip first data
                ERAsub.ddd = ERAsub.ddd(5:end);
                ERAsub.temp = ERAsub.temp(5:end,:,:);
                ERAsub.rh = ERAsub.rh(5:end,:,:);
                ERAgridTvP.era.time = ERAgridTvP.era.time(5:end);
                ERAgridTvP.era.T2 = ERAgridTvP.era.T2(5:end,:,:);
                ERAgridTvP.era.pvap = ERAgridTvP.era.pvap(5:end,:,:);
            end
            if(str2double(datestr(t,'mm'))==2 && leapyear(str2double(datestr(t,'yyyy')))==1)
                disp('29th Feb')
                ERAsub.ddd = [ERAsub.ddd; ERAgridrad.era.time(end-3:end)];
                ERAsub.temp = cat(1,ERAsub.temp, ERAsub.temp(end-3:end,:,:));
                ERAsub.rh = cat(1,ERAsub.rh, ERAsub.rh(end-3:end,:,:));
            end
        end
        if(sum((ERAsub.ddd)-(ERAgridrad.era.time)+(ERAgridTvP.era.time)-(ERAsub.ddd)) ~= 0)
            error('------- Still not syncron------')
        end
        tid = ERAgridrad.era.time;
        
        % reset vectors which are written to file
        SW = ones(nrows,ncols,length(tid),'single')*NaN;
        LW = ones(nrows,ncols,length(tid),'single')*NaN;
        
        %loop over sub-daily time step (6h)
        for h=1:length(tid)%h=27%(length(tid)-4):1:length(tid)%1:length(tid) ____________________________________________________________________________________________ sett tilbake _______________________________________________________________________-
            
            %disp([datestr(tid(h))])
            %disp(['lagging time ',datestr(tid(h)-.25)])
            
            % load Solar-Geometry
            ddd=doy(tid(h)-.25);
            [~,~,~,hh,~,~]=datevec(tid(h)-.25);
            hhh = num2str(100+hh); hhh=hhh(2:3);
            SUN=load([sun_path,num2str(ddd,'%03i'),'_hh',hhh,'.mat']);
            
            % interpolate ERA to Sval1000m
            SWgrid = interpERA2SVAL_simple(eraLON,eraLAT,ERAgridrad.era.SWsfc_6h(h,:,:),long,lat);
            SWtoa = interpERA2SVAL_simple(eraLON,eraLAT,ERAgridrad.era.SWtoa_6h(h,:,:),long,lat);
            SWtoa(SWtoa<.02)=0; % some noise occur ??
            
            % compute m
            m = 1./cosd(SUN.zen);
            
            % compute kt
            kt = SWgrid./SWtoa;
            kt(kt>1)=0;  % set to zero when SWgrid is higher, not good... inf occurs too.
            kt(isfinite(kt)==0)=0; % and nans
            
            % compute kd (diffuse fraction: after Ruiz-Arias et al 2010 eq. (16)
            kd = 0.952 - 1.041.*exp(-exp((2.300-4.702.*kt)));
            
            % partion in diffuse and direct
            SWeradif = kd.*SWgrid;
            SWeradir = SWgrid-SWeradif;
            
            % correct SWdir for elevation change
            SWsubdir = SWeradir;
            Z_correction = kt>=0.65; % ONLY CLEAR SKY, test if kt>0.65 Perez et al 1990
            k = -log(SWeradir./SWtoa)./m; % extintion coeff % after Ruiz-Arias et al 2010 eq. (1)
            SWsubdir(Z_correction) = SWeradir(Z_correction).*exp(-(k(Z_correction).*dP(Z_correction))./(P0*cosd(SUN.zen(Z_correction)))); % Kumar 1997 eq 8, 9
            %SWdir_dz = SWdir + SWdir*(1-exp(-k.*dz.*cosd(Z))); % after Ruiz-Arias et al 2010 eq. (2)
            %figure,imagesc(SWdir_dz-SWdir),title('altitude correction to SW_{dir} (SWdir_dz-SWdir)'),colorbar,axis xy
            S1=SWsubdir;
            % incidence angle on sloped surface
            %cosd_Ia = cosd(Z).*cosd(slope) + sind(Z).*cos(A-aspect); % after eq.(7) Fiddes & Gruber (2014)
            
            % Self sahdowing eq (9) in Fiddes & Gruber (2014)
            SWsubdir = SWsubdir .* (SUN.ia./cosd(SUN.zen));
            SWsubdir(SUN.zen>=89.5)=0; %S2=SWsubdir;
            
            % cast sahdow (binary mask)
            SWsubdir = SWsubdir .* SUN.cs;
            
            % Diffusive part  % after eq.(10) Fiddes & Gruber (2014)
            SWsubdif = SWeradif .* Vf;
            
            % final SW in terrain
            SWsub = SWsubdir + SWsubdif;
            SWsub(isnan(SWsub))=0;
            SWsub(SWsub<0)=0;
            
            
            % LONG WAVE RADATION
            % interpoate
            LW_grid = interpERA2SVAL_simple(eraLON,eraLAT,ERAgridrad.era.LWsfc_6h(h,:,:),long,lat);
            T2_grid = interpERA2SVAL_simple(eraLON,eraLAT,ERAgridTvP.era.T2(h,:,:),long,lat);
            pV_grid = 100*interpERA2SVAL_simple(eraLON,eraLAT,ERAgridTvP.era.pvap(h,:,:),long,lat); % convert from hPa to Pa
            
            % water vapor pressure
            %pV_grid = .01*RH_grid .* eso .* exp( (Lv/Rv) * ((1/To) - (1./T_grid(h,:,:)) ) );
            pV_sub = squeeze(ERAsub.rh(h,:,:) .* es0 .* exp( (Lv/Rv) * ((1/To) - (1./(ERAsub.temp(h,:,:))))));
            
            % clear sky emmisivity
            e_cl_grid = .23 + aa*(pV_grid./T2_grid).^(1/bb);
            e_cl_sub = squeeze(.23 + aa*(pV_sub./squeeze(ERAsub.temp(h,:,:))).^(1/bb));
            
            % all-sky emissivity
            e_as_grid = LW_grid ./ (sigma*(T2_grid.^4));
            e_diff = e_as_grid - e_cl_grid; % cloud component
            
            LW_sub = (e_cl_sub + e_diff) .* sigma .* (squeeze(ERAsub.temp(h,:,:)).^4);
            
            
            % store in monthly struct
            LW(:,:,h)=squeeze(LW_sub);
            SW(:,:,h)=squeeze(SWsub);
        end
        
        % test
        BADs =numel(SW)-(sum(sum(sum(isfinite(SW)))));
        if(BADs>0)
            disp([datestr(t),'  non-finite  ',num2str(BADs)])
        end
        
        %SAVE monthly files
        %file = [out_path,datestr(tid(1),'yyyymm'),'.mat'];
        %save(file,'tid','LW','SW','infoRAD');
        write_RADsub_file(out_path,tid,SW,LW,infoRAD)
    else
        disp([out_path,datestr(t,'yyyymm'),' exists, omitting...'])
    end %exits file
end

%matlabpool close
% delete(poolobj);
    
    
    
