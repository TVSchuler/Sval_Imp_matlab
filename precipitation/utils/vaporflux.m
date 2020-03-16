function WV = vaporflux(rho,q,u,v,Hw)
% function to calculate the column water vapor flux
% TVS, June 2014
    WV = rho*q*sqrt(u.^2+v.^2)*Hw; 
end