function res = linear_transition(x,xm,sigma)
    res                 = zeros(size(x));
    res(x>xm+sigma/2)   = 1;
    idx                 = find(x>xm-sigma/2&x<=xm+sigma/2);
    res(idx)            = (x(idx)-xm+sigma/2)/sigma;
end