function [lat,lon,p] = getERAtmp(tmpname,area)
    % get data
    tmp         = load(tmpname);

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
    p           = double(tmp.era.levelist);
end