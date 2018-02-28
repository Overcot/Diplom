function [ x ] = Boundary2(x0, u)
    global Ds Dt s t alpha gamma SSB_0 k L1;
    x(s, 1) = x0;
    for time=t(1:end - 1)
        SSB = trapz(gamma(1:end).*x(1:end, time))*Ds
        x(1, time) = L1/(1+exp(-k * (SSB - SSB_0)))
        x(s(2:end), time + 1) = x(s(2:end), time) + Dt*((alpha(s(2:end)-1)-1)'.*x(s(2:end)-1, time) - u(s(2:end)-1, time)) - (Dt/Ds)*(x(s(2:end), time)-x(s(2:end)-1, time));
    end
    x(1, end) = trapz(gamma(1:end).*x(1:end, t(end)))*Ds;
    
end