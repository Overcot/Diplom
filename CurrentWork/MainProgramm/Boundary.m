function [ x ] = Boundary(x0, u)
    global Ds Dt s t 
    global gamma mu;
    global a b allee;
    x(s, 1) = x0;
    for time=t(1:end - 1)

        %x(1, time) = recruitmentFunction('Anna', time, x, [Ds, Dt]);
        x(1, time) = recruitmentFunction('Tahvonen', time, x, [Ds, Dt]);

        % Potapov
        for class = s(1:end-1)
            %(x(class+1, time+1) - x(class, time))/Dt = -(mu(class+1)+mu(class))/(2*Dt) * (x(class+1, time+1)+x(class, time))/(2*Dt) - (u(s+h,t+h)+u(s,t))/(2*Dt) * (x(class+1, time+1)+x(class, time))/(2*Dt);
            upper = ( (mu(class + 1) + mu(class))/(2*Dt) * x(class, time)/(2*Dt) + (u(class + 1, time + 1) + u(class, time))/(2*Dt) * x(class, time)/(2*Dt) - x(class, time)/Dt);
            lower = (1/Dt + (mu(class + 1) + mu(class))/(2*Dt) + (u(class + 1, time + 1) + u(class, time))/(2*Dt));
            x(class+1, time+1) = - upper/lower;
        end
    end
    %x(1, t(end)) = recruitmentFunction('Anna', t(end), x, [Ds, Dt]);
    x(1, t(end)) = recruitmentFunction('Tahvonen', time, x, [Ds, Dt]);
    
end
