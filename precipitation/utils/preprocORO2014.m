function [oro,LSoro] = preprocORO2014(oro_file,g,area,region,X,Y)
    tmp             = load(oro_file);
    LSoro           = tmp.nc.z.data./g;
    LSlon           = tmp.nc.longitude.data;
    LSlat           = tmp.nc.latitude.data;

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
    switch region
        case 'svalbard'
            % streches over several UTMzones: we fix it to a single
            % zone: UTM33X --> use modified function
            [LSX,LSY,~]     = deg2utm33(LSLAT(:),LSLON(:));   
        otherwise
            % otherwise use this one
            [LSX,LSY,~]     = deg2utm(LSLAT(:),LSLON(:));
    end

    LSX1            = reshape(LSX,length(LSlat),length(LSlon));
    LSY1            = reshape(LSY,length(LSlat),length(LSlon));
    oro             = griddata(LSX1,LSY1,LSoro,X,Y,'cubic');
    oro(oro<0)      = 0; % no depressions here
    % also the LSoro as output, clip to the area
    LSoro           = LSoro(2:end-1,2:end-1);
 
end