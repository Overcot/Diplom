function [ x ] = Boundary(x0, u)
    global Ds Dt s t 
    global gamma mu;
    global a b allee;
    x(s, 1) = x0;
    for time=t(1:end - 1)
        %{
        prev = 1;
        while prev~=x(1,time)
            prev = x(1,time);
            x(1, time) = recruitmentFunction('Anna', time, x, [Ds, Dt]);
        end
        %}
        x(1, time) = recruitmentFunction('Anna', time, x, [Ds, Dt]);
        %x(s(2:end), time + 1) = x(s(2:end), time) + Dt*((alpha(s(2:end)-1)-1)'.*x(s(2:end)-1, time) - u(s(2:end)-1, time)) - (Dt/Ds)*(x(s(2:end), time)-x(s(2:end)-1, time));

        % Potapov
        for class = s(1:end-1)
            %(x(class+1, time+1) - x(class, time))/Dt = -mu(class) * (x(class+1, time+1)+x(class, time))/2;
            %x(class+1,time+1)/Dt + mu(class)*x(class+1, time+1)/2 = x(class, time)/Dt - mu(class)*x(class, time)/2;
            x(class+1, time+1) = (x(class, time)/Dt - mu(class)*x(class, time)/2)/(1/Dt + mu(class)/2);
            
        end
    end
    x(1, t(end)) = recruitmentFunction('Anna', t(end), x, [Ds, Dt]);
    
end
