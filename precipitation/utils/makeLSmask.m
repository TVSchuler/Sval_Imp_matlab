function LSmask = makeLSmask(oro_file,g,area,region)
    tmp             = load(oro_file);
    LSoro           = tmp.nc.z.data./g;
    LSlon           = tmp.nc.longitude.data;
    LSlat           = tmp.nc.latitude.data;

    inlat           = find(LSlat>=area(3)&LSlat<=area(4));
    inlat           = inlat(1)-3:inlat(end)+3;

    inlon           = find(LSlon>=area(1)&LSlon<=area(2));
    inlon           = inlon(1)-8:inlon(end)+6;

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
    
    % project onto a regular grid
    minX=min(min(LSX1));
    maxX=max(LSX1(:));
    minY=min(LSY1(:));
    maxY=max(LSY1(:));
    
    % regular 20km spacing in x and y
    x=minX:2e4:maxX;
    y=minY:2e4:maxY;
    
    [X,Y]=meshgrid(x,y);
    
    oro             = griddata(LSX1,LSY1,LSoro,X,Y,'cubic');
    % distance mapfrom the 200m contour...
    d200 = bwdist(oro>=200);
    % mask= 1 everywhere up to 200km distance from the 200m contour
    % -->dist = 10 (*20km spacing)
    mask = d200>0&d200<=10;
    
    % reproject the mask onto the LSoro grid used...
%     LSlon = LSlon(9:end-6);
%     LSlat = LSlat(4:end-3);
LSX2 = LSX1(4:end-3,9:end-6);
LSY2 = LSY1(4:end-3,9:end-6);

    
    LSmask         = griddata(X,Y,double(mask),LSX2,LSY2,'nearest');
    LSmask = floor(LSmask);
end
    
    