function w           = mixingratio_weighting(q,z)
% calculates weights w for vertical averaging over elevation z based on moisture content q
    [~,nlev,~,~]     = size(z);
    dz               = diff(z,1,2);
    % trapezoidal integration of q
    qdz(:,1,:,:)         = q(:,1,:,:).*dz(:,1,:,:)/2;
    for i=2:nlev-1
        qdz(:,i,:,:)     = q(:,i,:,:).*(dz(:,i,:,:)+dz(:,i-1,:,:))/2;
    end
    qdz(:,nlev,:,:)      = q(:,nlev,:,:).*dz(:,nlev-1,:,:)/2;
    q_tot            = sum(qdz,2);
    % preallocate w
    w                    = z.*0;
    w(:,1,:,:)           = qdz(:,1,:,:)./q_tot;
    for i=2:nlev-1
        w(:,i,:,:)       = qdz(:,i,:,:)./q_tot;
    end
    w(:,nlev,:,:)        = qdz(:,nlev,:,:)./q_tot;
end