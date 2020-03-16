function precip  = apply_LTmodel2014(A,hprime,ncx,ncy,ncol,nrow,cellsize)

    % expand the A structure...
    sens            = A.sens;
    rho             = A.rho;
    q               = A.q;
    hw              = A.hw;
    u               = A.U;
    v               = A.V;
    Nmsq            = A.Nmsq;
    tauf            = A.tauf;
    tauc            = A.tauc;

    % wave numbers

    % wavenumber increment
    dkx           = 2*pi/(ncx*cellsize);
    dky           = 2*pi/(ncy*cellsize);
    % construct wavenumber arrays
    ncxh          = floor(ncx/2)+1;
    kx_1          = dkx*(0:ncxh-1);
    if(rem(ncx,2)==0)
        kx_2          = -fliplr(kx_1(2:end-1));
    else
        kx_2          = -fliplr(kx_1(2:end));
    end
    kx            = [kx_1,kx_2];

    ncyh          = floor(ncy/2)+1;
    ky_1          = dky*(0:ncyh-1);
    if(rem(ncy,2)==0)
        ky_2          = -fliplr(ky_1(2:end-1));
    else
        ky_2          = -fliplr(ky_1(2:end));
    end
    ky            = [ky_1,ky_2];

    [KX,KY]       = meshgrid(kx,ky);
    
    % this is a bit ugly... we have to turn around the v-component of
    % wind due to the sign convention: v is positive for a wind from S to
    % N, but the spectral method advects precip along the matrix, i.e. for
    % a positive v that would mean from N to S...
    v             = -v;
    % continue...
    sigma         = u*KX+v*KY;
    % size(sigma)
    denom         = sigma.^2;
    % denom(abs(denom)<1e-18)=sign(denom(abs(denom)<1e-18))*1e-18;   % avoid division by zero!!
    denom(abs(denom)<1e-18)=1e-18;   % avoid division by zero!!
    mtermsq       = ((Nmsq-sigma.^2)./denom).*(KX.^2+KY.^2);

    mterm         = zeros(size(mtermsq));

    % propagating case
    pos           = find(mtermsq>=0);
    mterm(pos)    = complex(sign(sigma(pos)).*sqrt(mtermsq(pos)),zeros(size(pos)));
    % evanescent case:
    neg           = find(mtermsq<0);
    mterm(neg)    = complex(zeros(size(neg)),sqrt(-mtermsq(neg)));

    % consider delay time
    fctr          = 1./(...
        (1-1i*mterm*hw)...
        .*(1+1i*sigma*tauc)...
        .*(1+1i*sigma*tauf)...
        );

    % calculate precip in frequency space: eq 49 in S&B04
    p_cplx       = sens*rho*q.*hprime*1i.*sigma.*fctr;

    % bring it back
    ppp           = ifft2(p_cplx);
        % size(ppp)
    ppp           = ppp(1:nrow,1:ncol);
       
    precip        = real(ppp); % in [kg m^-2 s^-1]

% avoid negative values
% precip(precip<0)=0; % actually, we should preserve negatives in the
% orographic correction, but remove them from the final precip in the main
% program
end