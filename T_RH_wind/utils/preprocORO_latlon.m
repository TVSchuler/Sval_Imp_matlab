function oro = preprocORO_latlon(oro_file,g,area,lat,lon)
    tmp             = load(oro_file);
    LSoro           = tmp.nc.Geopotential.data./g;
    LSlon           = tmp.nc.lon.data;
    LSlat           = tmp.nc.lat.data;

    inlat           = find(LSlat>=area(3)&LSlat<=area(4));
    inlat           = inlat(1)-1:inlat(end)+1;

    inlon           = find(LSlon>=area(1)&LSlon<=area(2));
    inlon           = inlon(1)-1:inlon(end)+1;

    % cut out the ROI (add 1 value to each side...
    LSlon           = LSlon(inlon);
    LSlat           = LSlat(inlat);

    LSoro           = LSoro(inlat,inlon);
    % mesh and convert
    [LSLON,LSLAT]   = meshgrid(LSlon,LSlat);
    [LSX,LSY,~]     = deg2utm33(LSLAT(:),LSLON(:));
    [LON,LAT]       = meshgrid(lon,lat);
    [X,Y,~]         = deg2utm33(LAT(:),LON(:));

%     LSX1            = reshape(LSX,length(LSlat),length(LSlon));
%     LSY1            = reshape(LSY,length(LSlat),length(LSlon));
    F               = scatteredInterpolant(LSX,LSY,LSoro(:));
    oro             = F(X,Y);
    oro(oro<0)      = 0; % no depressions here
    oro             = reshape(oro,length(lat),length(lon));
    
end