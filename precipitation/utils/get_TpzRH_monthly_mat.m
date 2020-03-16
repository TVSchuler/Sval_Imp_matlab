function [T,p,z,RH,U,V,time,lon,lat] = get_TpzRH_monthly_mat(plfile,area,g)

    % get met data
    tmp         = load(plfile);

    % date vector for this file
    time        = tmp.era.time;

    % get lat and lon
    lon         = tmp.era.longitude;
    lat         = tmp.era.latitude;

    % extract the region defined by area
    if(isempty(area))   % default is entire area...
        id1         = 1:numel(lon);
        id2         = 1:numel(lat);
    else
        id1         = find(lon>=area(1)&lon<=area(2));  % longitudes go from -180 to 180
        id2         = find(lat>=area(3)&lat<=area(4));  % latitude goes from 90 to -90
    end
    lon         = lon(id1);
    lat         = lat(id2);
    % use the idx to extract only the region corresponding to area
    % the dimensions are: ens, time, isobaric, lat, lon !!!

    % get geopotential height and temperatures, RH and wind components
    z           = tmp.era.Z(:,:,id2,id1)./g;
    T           = tmp.era.T(:,:,id2,id1);
    RH          = tmp.era.RH(:,:,id2,id1);
    U           = tmp.era.U(:,:,id2,id1);
    V           = tmp.era.V(:,:,id2,id1);

    % make up p from constant info
    % p           = repmat([1000,925,850,700],[length(idx),numel(lat),numel(lon)]);
    p           = double(tmp.era.levelist);
end
