function [x, x0, u, Ds, Dt, S_steps, T_steps, s, t, gamma, mu] = startInit(precision, L, T, gammaData, muData, x0Data, uData)
    Ds = precision;
    Dt = precision;
    S_steps= L/Ds;
    T_steps = T/Dt;
    s=1:S_steps+1;
    t=1:T_steps+1;
    
    gamma = gammaData;
    mu = muData;
    
    x = zeros(size(s, 2), size(t, 2));
    x0(s(1:(S_steps/L):end), 1) = x0Data;
    % x(s(1:end), 1) = interp1(s(1:(S_steps/L):end),x0(s(1:(S_steps/L):end)),s(1:end));
    
    u = zeros(size(s,2),size(t,2));
    u = uData;
end

