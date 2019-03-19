function [ x ] = Boundary(x0, u)
    global Ds Dt s t 
    global gamma mu;
    x(s, 1) = x0;
    for time=t(1:end - 1)

        x(1, time + 1) = recruitmentFunction('Anna', time, x, [Ds, Dt]);
        % не очень правильное тановен
        %x(1, time) = recruitmentFunction('Tahvonen', time, x, [Ds, Dt]);

        % Potapov
        for class = s(1:end-1)
            %(x(class+1, time+1) - x(class, time))/Dt = -(mu(class) * x(class, time) - u(class, time) * x(class, time));
            x(class+1, time+1) = Dt*(-mu(class) * x(class, time) - u(class, time) * x(class, time) + x(class, time)/Dt);
        end
    end
    %x(1, t(end)) = recruitmentFunction('Anna', t(end), x, [Ds, Dt]);
    %x(1, t(end)) = recruitmentFunction('Tahvonen', t(end), x, [Ds, Dt]);
    
end
