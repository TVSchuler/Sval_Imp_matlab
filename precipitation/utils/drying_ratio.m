function [DR]=drying_ratio(P,WV,u,v,cellsize)
    % function to calculate the drying ratio according to Smith&Evans (2007),
    % Appendix C; the integral in equation C5 is calculated in 2D for 8
    % wind-direction classes
    % TVS, June 2014
    
    % do not allow negative P
    P               = max(P,0);
    
    % determine size of matrix
    [nrow,ncol]     = size(P);
    % number of elements    
    nP              = numel(P);
    
    % wind direction
    rad2deg         = 180/pi;
    wdir            = atan2(u,v)*rad2deg; % the direction towards which the wind blows...
    wdir            = mod(wdir+22.5,360);   
    % output from atan2d is in degrees [-180, 180], convert to [0 360]
    % and turn it by 22.5 deg such that limits of classes are at 
    % 22.5, 67.5...and not at 45, 90 etc
    % atan2d was introduced in R2012b, to keep compatibility with older
    % versions, we use atan2 and do the rad2deg conversion

    % classify into 8 classes:
    wclass          = floor(wdir/45);
    
    % sum up precip in direction of wind, 8 classes
    switch wclass
        case 0 % wind from S
            DR = flipud(cumsum(flipud(P)))/WV*cellsize; % P*cellsize^2/WV*cellsize
        case 1 % wind from SW
            A  = flipud(P);
            % initialise DR
            DR = zeros(size(A))*NaN;
            for d =-nrow+1:ncol-1   % for all diagonals, 
                % we start from lower left and successively overwrite the
                % erroneously addressed elements in the upper corner
                
                % for negative d
                if d<=0
                    % index to the diagonale
                    idx     = (-d+1:nrow+1:nP);
                    DR(idx) = cumsum(A(idx))/WV/sqrt(2)*cellsize;
                else % for positive d
                    idx     = (d*nrow+1:nrow+1:nP);
                    DR(idx) = cumsum(A(idx))/WV/sqrt(2)*cellsize;
                end
            end
            % flip it back...
            DR = flipud(DR);
        case 2 % wind from W
            DR = cumsum(P,2)/WV*cellsize;
        case 3 % wind from NW
            DR = zeros(size(P))*NaN;
            for d =-nrow+1:ncol-1   % for all diagonals
                % for negative d
                if d<=0
                    % index to the diagonale
                    idx     = (-d+1:nrow+1:nP);
                    DR(idx) = cumsum(P(idx))/WV/sqrt(2)*cellsize;
                else % for positive d
                    idx     = (d*nrow+1:nrow+1:nP);
                    DR(idx) = cumsum(P(idx))/WV/sqrt(2)*cellsize;
                end
            end
        case 4 % wind from N
            DR = cumsum(P)/WV*cellsize;
        case 5 % wind from NE
            A = fliplr(P);
            % initialise DR
            DR = zeros(size(A))*NaN;
            for d =-nrow+1:ncol-1   % for all diagonals
                % for negative d
                if d<=0
                    % index to the diagonale
                    idx     = (-d+1:nrow+1:nP);
                    DR(idx) = cumsum(A(idx))/WV/sqrt(2)*cellsize;
                else % for positive d
                    idx     = (d*nrow+1:nrow+1:nP);
                    DR(idx) = cumsum(A(idx))/WV/sqrt(2)*cellsize;
                end
            end
            % flip it back...
            DR = fliplr(DR);
        case 6 % wind from E
            DR = fliplr(cumsum(fliplr(P),2))/WV*cellsize;
        case 7 % wind from SE
            A = flip(flip(P,2));
            DR = zeros(size(A))*NaN;
            for d =-nrow+1:ncol-1   % for all diagonals
                % for negative d
                if d<=0
                    % index to the diagonale
                    idx     = (-d+1:nrow+1:nP);
                    DR(idx) = cumsum(A(idx))/WV/sqrt(2)*cellsize;
                else % for positive d
                    idx     = (d*nrow+1:nrow+1:nP);
                    DR(idx) = cumsum(A(idx))/WV/sqrt(2)*cellsize;
                end
            end
            % flip it back:
            DR = flip(flip(DR,2));
        otherwise
            error('something went wrong with wind classes in drying_ratio.m')
    end

end