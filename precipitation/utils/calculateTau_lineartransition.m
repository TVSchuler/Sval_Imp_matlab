function [tauc,tauf]    = calculateTau_lineartransition(hw_m,Tmean)
    % to calculate tau, we use hw and our v_f parameterization
    vmin          = 1.; % m/s fall speed of snow
    vmax          = 1.5; % m/s fall speed of rain

    % v_f from linear transition
    xm            = 273; % K centre of transition
    sigma         = 5;  % K width of transition
    vscale        = linear_transition(Tmean,xm,sigma);
    v_f           = (vmax-vmin)*vscale + vmin;

    tauf          = hw_m./v_f;
    tauf          = max(tauf,10);

    tauc          = tauf;
end