function [zzz,X,Y,cellsize,nrow,ncol] = preprocDEM2014(demname)
    [zzz, ncol, nrow, xll, yll, cellsize, ~] = asciigridread(demname);
    zzz(zzz<0)  = 0; % we do not have depressions
    % zzz         = flipud(zzz); why should we do that??
    xkoo        = xll+(0:ncol-1)*cellsize;
    ykoo        = yll+(0:nrow-1)*cellsize;
    % flip y around: 1rst pnt is upper left corner
    ykoo = flipud(ykoo');

    [X,Y]       = meshgrid(xkoo,ykoo);

end