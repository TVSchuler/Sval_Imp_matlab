function Z = interpERA2SVAL(eraLON,eraLAT,eraZ,x,y)
% interpolate ERA from long lat grid to Svalbard 1000m UTM33X wgs84
% addpath('/uio/kant/geo-natg-u1/torbjoos/mToolBox/','/uio/kant/geo-natg-u1/torbjoos/mToolBox/m_map14/m_map/');


%% mesh
[long,lat] = meshgrid(double(eraLON),double(eraLAT));
[X,Y] = meshgrid(x,y);

%% convert too utm
longlim = [5 37];
latlim  = [74 83];
m_proj('UTM','ellipsoid','wgs84','longitude',longlim,'latitude',latlim,'zone',33);
% convert
[eraX,eraY] = m_ll2xy(long,lat);

%% interpolate
F = scatteredInterpolant(eraX(:),eraY(:),double(eraZ(:)));
Z = F(X(:),Y(:));

Z=reshape(Z,size(X));

