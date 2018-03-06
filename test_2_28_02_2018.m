clear; clc;
precision = 0.1;
L = 51;
T = 250;
Ds = precision;
Dt = precision;
S_steps= L/Ds;
T_steps = T/Dt;
s=1:S_steps+1;
t=1:T_steps+1;
gamma(1:S_steps/2) = 0;
gamma(S_steps/2+1:s(end)) = 2/L;
x(s(1:end), 1) = 1;
mu(s(1:end)) = (s-1)/L*Ds;
u(s(1:end),t(1:end)) = 0;

for time=t(1:end)
    u(s(1:end),time) = -1/((time+1)*Dt)-(s-1)/L*Ds;
end

for time=t(1:end - 1)
    x(s(2:end-1), time + 1) = x(s(2:end-1), time) - (Dt/(2*Ds))*(x(s(2:end-1)+1, time) - x(s(2:end-1)-1, time)) - Dt*(mu(s(2:end-1))'.*x(s(2:end-1), time) + u(s(2:end-1),time).*x(s(2:end-1), time));
    x(s(end), time + 1) = x(s(end), time) - (Dt/Ds)*(x(s(end), time) - x(s(end)-1, time)) - Dt*(mu(s(end))'.*x(s(end), time) + u(s(end),time).*x(s(end), time));
    x(1,time+1) = trapz(gamma(1:end)'.*x(1:end, time)*Ds);
end

plot(x(1, t))
