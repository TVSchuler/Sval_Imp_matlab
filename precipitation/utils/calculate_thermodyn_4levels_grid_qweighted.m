function [mean_gamma, mean_Gamma_m, Nmsqu, Hw, senscoef, qbase, mean_T, w1]     = calculate_thermodyn_4levels_grid_qweighted(T,z,p,RH,g)
    % Gamma_m (following Durran & Klemp 1982, eq 36)
    T_0                = 273.15;                          % 0C in K
    c_p                = 1004;                            % heat capacity of air [J kg-1 K-1]
    % c_p                = 4181;                            % heat capacity of water [J kg-1 K-1]
    R                  = 287;                             % gas constant for dry air [J kg-1 K-1]
    Rv                 = 461;                             % Rv = gas constant for vapor [J kg-1 K-1]
    eps                = 0.622;                           % ratio R/Rv [-]; 
    p_0                = 1000;                            % [hPa]
    e_s0               = 6.112;                           % e_s at T=0 C [hPa];

    Gamma_d            = g/c_p;                           % dry adiabatic lapse rate [K m-1]
    
    % temperature gradient
    dTdz               = diff(T,1,2)./diff(z,1,2);

% the following expressions are from Stone&Carlson 1979
% L as a function of T
    L_T                = (2510 - 2.38*(T-T_0))*1e3;    % [J kg-1]
    L_T                = squeeze(L_T);
    % at intermediate levels
%     L_T_m               = [(L_T(:,1)+L_T(:,2))/2,(L_T(:,2)+L_T(:,3))/2];

    % calculate e_s from T: Bolton 1980
    e_1                = e_s0*exp(17.67*(T-T_0)./(T-29.65));
    % actual vapor pressure, RH is in [%]...
    e_2                = RH/100.*e_1;
    
    q                  = zeros(size(T));
    qs                 = zeros(size(T));
    Gamma_m            = zeros(size(T));
    
    for i=1:numel(p)
        qs(:,i,:,:)    = eps*e_1(:,i,:,:)./(p(i)-e_1(:,i,:,:));                  % saturation mixing ratio
        q(:,i,:,:)     = eps*e_2(:,i,:,:)./(p(i)-e_2(:,i,:,:));
    end

    gamma              = -dTdz;

    for i=1:numel(p)
        % --> AMS glossary: moist adiabatic lapse rate
        Gamma_m(:,i,:,:) = g.*(1+(L_T(:,i,:,:).*qs(:,i,:,:))...
            ./(R.*T(:,i,:,:)))./(c_p+(L_T(:,i,:,:).^2.*qs(:,i,:,:).*eps)...
            ./(R.*T(:,i,:,:).^2));
    end

    % use the mixing ratio for weighting
    w1                 = mixingratio_weighting(q,z);
    
    w2                 = inter_mixingratio_weighting(q,z);
% not used since we have modified calculation of mean_gamma

        mean_gamma    = sum(gamma.*w2,2);  % the q-weighted version is
    %     abandoned, it is the max lapse rate that represents the lowest
    %     stability --> governs overall stability
    % we do not consider the 1000hPa which is typically influenced by the
    % surface, but we want to characterise the air mass as a whole
%     mean_gamma         = max(gamma(:,2:3,:,:),[],2); 
%     mean_gamma         = gamma(:,2,:,:);
%     mean_gamma         = max(mean_gamma,1e-3);
%     mean_gamma         = min(mean_gamma,7e-3);
mean_gamma               = -(T(:,4,:,:)-T(:,2,:,:))./(z(:,4,:,:)-z(:,2,:,:));
    
%     mean_Gamma_m       = sum(Gamma_m.*w1,2);
    mean_Gamma_m       = Gamma_m(:,4,:,:)-Gamma_m(:,2,:,:)./(z(:,4,:,:)-z(:,2,:,:));
    mean_T             = sum(T.*w1,2);
    mean_L_T           = sum(L_T.*w1,2);
%     mean_Gamma_m       = Gamma_m(:,2,:,:);
    mean_T             = T(:,1,:,:);
    mean_L_T           = L_T(:,1,:,:);
%     mean_Gamma_m       = mean(Gamma_m(:,2:3,:,:),2);
%     mean_Gamma_m       = max(mean_Gamma_m,5e-3);
%     mean_Gamma_m       = min(mean_Gamma_m,9e-3);
    
%     mean_T             = mean(T(:,2:3,:,:),2);
%     mean_L_T           = mean(L_T(:,2:3,:,:),2);
    
    Nmsqu              = g./mean_T.*(mean_Gamma_m - mean_gamma); % = equ3  which is an OK approximation (D&K1982)!
    
    % if Nmsqu<0 --> unstable stratification --> we set Nmsqu=0 and calculate
    % hw and sens correspondingly    
    Nmsqu              = max(Nmsqu,0);
    % set an upper limit for Nmsqu
    Nmsqu              = min(Nmsqu,1e-4);
    % recalculate mean_gamma, consistent with new Nmsqu
    mean_gamma         = -Nmsqu.*mean_T/g+mean_Gamma_m;
    if(min(mean_gamma(:))<=0)
        disp('warning: negativ lapse rate!')
    end
    Hw                 = Rv*T(:,1,:,:).^2./(mean_L_T.*mean_gamma); % minus in S&B04 ??? --> definition of gamma!
%     Hw                 = Rv*mean_T.^2./(mean_L_T.*mean_gamma); % minus in S&B04 ??? --> definition of gamma!
    senscoef           = mean_Gamma_m./mean_gamma;  
    % only lowermost value for qs
%     qbase              = mean(qs(:,2:3,:,:),2); % saturation mixing ratio
    qbase              = qs(:,1,:,:); % saturation mixing ratio
% qbase = sum(qs.*w1,2);
end