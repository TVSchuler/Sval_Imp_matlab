function [zzz,X,Y,cellsize,nrow,ncol,FX,FY] = preprocDEM(demname)
    [zzz, ncol, nrow, xll, yll, cellsize, ~] = asciigridread(demname);
    zzz(zzz<0)  = 0; % we do not have depressions
    % zzz         = flipud(zzz); why should we do that??
    xkoo        = xll+(0:ncol-1)*cellsize;
    ykoo        = yll+(0:nrow-1)*cellsize;
    % flip y around: 1rst pnt is upper left corner
    ykoo = flipud(ykoo');

    [X,Y]       = meshgrid(xkoo,ykoo);

    % slope of the terrain 
    [FX,FY]     = gradient(zzz,cellsize,cellsize);
    % % gradient is applied on the matrix, 
    % % so positive val means increasing z with rownr, 
    % % but we would like to have pos. val for increase with North-coordinate
    FY          = -FY; 

    % set boundaries of slope to zero
    FY(1,:)     = 0;
    FY(end,:)   =  0;

    FX(:,1)     = 0;
    FX(:,end)   = 0;
end