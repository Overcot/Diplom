clear; clc;
precision = 1;
L = 5;
T =1 ;
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
    u(s(1:end),time) = -1/(time)*Dt - (s-1)/L*Ds;
end
%{
A = -1/(2*Ds);
B = 1/Dt;
C = 1/(2*Ds);
%}

for time=t(1:end - 1)
    alpha(s) = zeros(size(s,2),1);
    beta(s) = zeros(size(s,2),1);
    
    %x(1,time+1) = gamma(1:end)*x(1:end, time)*Ds;

    %считаем самую правую точку с помощью €вного левого уголка
    %x(s(end), time + 1) = x(s(end), time) - (Dt/Ds)*(x(s(end), time) - x(s(end)-1, time)) - Dt*(mu(s(end))'.*x(s(end), time) + u(s(end),time).*x(s(end), time));
    
    %Ќе€вна€ схема(метод прогонки)
    %{
    for class=s(2:end-1)
        Fij = -mu(class)*x(class, time)-u(class, time)*x(class, time) + x(class, time)/Dt;
        if class == s(end)-1
            beta(class) = (Fij - A*beta(class-1))/(A*alpha(class-1)+B);
        elseif class == 2
            alpha(class) = -C/B;
            beta(class) = Fij/B;
        else
            alpha(class) = -C/(A*alpha(class-1)+B);
            beta(class) = (Fij - A*beta(class-1))/(A*alpha(class-1)+B);
        end
    end
    %}
    %{
    temp = zeros(size(s,2)-2,size(s,2)-2);
    tempf = zeros(size(s,2)-2,1);
    A = -1/(2*Ds);
    C = 1/(2*Ds);
    for class=s(2:end-1)
        classx = class-1
        if class == 2
            temp(classx, classx) = 1/Dt + mu(class) + u(class, time + 1);
            temp(classx, classx + 1) = 1/(2*Ds);
            tempf(classx) = x(class, time)/Dt - A*x(1, time);
        elseif class == s(end-1)
            temp(classx, classx-1) = -1/(2*Ds);
            temp(classx, classx) = 1/Dt + mu(class) + u(class, time + 1);
            tempf(classx) = x(class, time)/Dt - C*x(s(end), time);
        else
            temp(classx, classx - 1) = -1/(2*Ds);
            temp(classx, classx) = 1/Dt + mu(class) + u(class, time + 1);
            temp(classx, classx + 1) = 1/(2*Ds);
            tempf(classx) = x(class, time)/Dt;
        end
    end
    temp;
    tempf;
    linsolve(temp, tempf)
    %}
    %{
    for class=s(2:end-1)
        A = -1/(2*Ds);
        B = 1/Dt + mu(class) + u(class, time + 1);
        C = 1/(2*Ds);
        Fij = x(class, time)/Dt;
        if class == s(end)-1
            beta(class) = (Fij - A*beta(class-1))/(A*alpha(class-1)+B);
        elseif class == 2
            alpha(class) = -C/B;
            beta(class) = (Fij-A*x(1, time))/B;
        else
            alpha(class) = -C/(A*alpha(class-1)+B);
            beta(class) = (Fij - A*beta(class-1))/(A*alpha(class-1)+B);
        end
    end
    %}
    
    for class=s(1:end)
        A = -1/(2*Ds);
        B = 1/Dt + mu(class) + u(class, time + 1);
        C = 1/(2*Ds);
        Fij = x(class, time)/Dt;
        if class == s(end)
            beta(class) = (Fij - A*beta(class-1))/(A*alpha(class-1)+B);
        elseif class == 1
            alpha(class) = -C/B;
            beta(class) = (Fij)/B;
        else
            alpha(class) = -C/(A*alpha(class-1)+B);
            beta(class) = (Fij - A*beta(class-1))/(A*alpha(class-1)+B);
        end
    end
    x(s(end), time + 1) = beta(s(end));
    for class=s(end-1:-1:1)
        x(class, time + 1) = beta(class) + alpha(class) * x(class + 1, time + 1)
    end

    
    %{
    x(s(2:end-1), time + 1) = x(s(2:end-1), time) - (Dt/(2*Ds))*(x(s(2:end-1)+1, time) - x(s(2:end-1)-1, time)) - Dt*(mu(s(2:end-1))'.*x(s(2:end-1), time) + u(s(2:end-1),time).*x(s(2:end-1), time));
    x(s(end), time + 1) = x(s(end), time) - (Dt/Ds)*(x(s(end), time) - x(s(end)-1, time)) - Dt*(mu(s(end))'.*x(s(end), time) + u(s(end),time).*x(s(end), time));
    %}
    
    %считаем самую левую точку из краевого услови€ в задаче
    
    
    
    
end
% x(1,t(end)) = gamma(1:end)*x(1:end, t(end))*Ds;
%plot(x(1, t))
