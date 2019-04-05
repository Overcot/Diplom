function [ x ] = Boundary(x0, u)
    global Ds Dt s t;
    global gamma mu;
    % стандартный формат для матриц x, u - LxT
    % Пример: в первом столбце распределение по рыбам в момент времени 0, по всем возрастам
    x(s, 1) = x0;
    for time=t(1:end - 1)

        x(1, time + 1) = recruitmentFunction('Anna', time, x);
        %x(1, time) = recruitmentFunction('Tahvonen', time, x);

        % Potapov
        for class = s(1:end-2)
            %(x(class+1, time+1) - x(class, time))/Dt = -(mu(class) * x(class, time) - u(class, time) * (1-mu(class))* x(class, time));
            x(class + 1, time + 1) = Dt*(-mu(class) * x(class, time) - u(class, time) * (1 - mu(class)) * x(class, time) + x(class, time)/Dt);
        end
        class = s(end) - 1;
        lastClass = s(end);
        x(s(end), time + 1) = Dt*(-mu(class) * x(class, time) - u(class, time) * (1 - mu(class)) * x(class, time) + x(class, time)/Dt) + Dt*(-mu(lastClass) * x(lastClass, time) - u(lastClass, time) * (1 - mu(lastClass)) * x(lastClass, time) + x(lastClass, time)/Dt);
    end
    for time=t(1:end)
        for class=s(1:end)
            if x(class, time) < 10e-3
                x(class, time) = 0;
            end
        end
    end
end
