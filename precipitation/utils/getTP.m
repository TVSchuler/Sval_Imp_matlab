function [TP,TPdate,TPX,TPY]=getTP(sffile)
    % load TP
%     disp('load TP')
    tmp           = load(sffile); % TP6h(t,lat,lon),date6h, lat, lon, units
    TP            = tmp.eraTP.TP; % dims: time,lat,lon
    TPdate        = tmp.eraTP.time;
    % map the TP to the DEM
    % the lat,lon of the TP-grid
    TPlat         = tmp.eraTP.latitude;
    TPlon         = tmp.eraTP.longitude;

    [TPLON,TPLAT] = meshgrid(TPlon,TPlat);

    % convert to UTM33
    [TPX,TPY,~]      = deg2utm33(TPLAT(:),TPLON(:));
    TPX              = reshape(TPX,length(TPlat),length(TPlon));
    TPY              = reshape(TPY,length(TPlat),length(TPlon));
end