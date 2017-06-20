function [ x ] = Boundary(x0, u )
    global Ds Dt s t alpha f gamma;
    x(s, 1) = x0;
    for time=t(1:end - 1)
        x(1, time) = trapz(gamma(1:end).*x(1:end, time))*Ds + f(time);
        x(s(2:end), time + 1) = x(s(2:end), time) + Dt*((alpha(s(2:end)-1)-1)'.*x(s(2:end)-1, time) - u(s(2:end)-1, time)) - (Dt/Ds)*(x(s(2:end), time)-x(s(2:end)-1, time));
    end
    x(1, end) = trapz(gamma(1:end).*x(1:end, t(end)))*Ds + f(t(end));

end

