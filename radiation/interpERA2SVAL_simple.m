function Z = interpERA2SVAL_simple(eraLON,eraLAT,eraZ,LON,LAT)
% interpolate ERA from long lat grid to Svalbard 1000m UTM33X wgs84

%% interpolate
F = scatteredInterpolant(eraLON(:),eraLAT(:),double(eraZ(:)));
Z = F(LON(:),LAT(:));

Z=reshape(Z,size(LON));

