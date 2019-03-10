function [x, x0_new, u, Ds, Dt, S_steps, T_steps, s, t, gamma, mu] = startInit(precision, L, T, gammaValue, muValue, x0)
    
    Ds = precision;
    Dt = precision;
    S_steps= L/Ds;
    T_steps = T/Dt;
    s=1:S_steps+1;
    t=1:T_steps+1;
    
    gamma = gammaValue;
    mu = muValue;
    
    x = zeros(size(s,2),size(t,2));
    x0_new(s(1:(S_steps/L):end)) = x0;
    x(s(1:end), 1) = interp1(s(1:(S_steps/L):end),x0_new(s(1:(S_steps/L):end)),s(1:end));
    
    u = zeros(size(s,2),size(t,2));
    
end

