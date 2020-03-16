function w           = inter_mixingratio_weighting(q,z)
% similar as in mixingratio_weighting.m but for values at intermediate
% z-levels
    [~,nlev,~,~]         = size(z);
    dz               = diff(z,1,2);
    % trapezoidal integration of q
    qdz              = (q(:,1:nlev-1,:,:)+q(:,2:nlev,:,:))/2.*dz(:,1:nlev-1,:,:);
    
    q_tot            = sum(qdz,2);
    % preallocate w
    w                = qdz.*0;
    for i=1:nlev-1
        w(:,i,:,:)       = qdz(:,i,:,:)./q_tot;
    end
    
end