function [x, L, T, Ds, Dt, S_steps, T_steps, s, t, gamma, mu, p] = startInit(precision, L, T, gammaValue, muValue, x0)
       
    Ds = precision;
    Dt = precision;
    S_steps= L/Ds;
    T_steps = T/Dt;
    s=1:S_steps+1;
    t=1:T_steps+1;
    
    gamma = gammaValue;
    mu = muValue;
    p(t(1:end)) = ((t-1)*Dt)./((t-1)*Dt+1);

    x = zeros(size(s,2),size(t,2));
    x(s(1:end), 1) = x0;
    
end

