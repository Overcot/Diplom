clc; clear;
% вправо время увеличивается
% вниз класс увеличивается
%% Params of Setki
global Ds Dt s t alpha f gamma umin umax;

L = 7;
T = 8;
Ds = 0.1;
Dt = 0.1;
S_steps= L/Ds;
T_steps = T/Dt;
s=[1:S_steps+1];
t=[1:T_steps+1];
alpha(s(1:(S_steps/L):end)) = [1 0.95 0.85 0.8 0.7 0.5 0.3 0];
alpha(s) = interp1(s(1:(S_steps/L):end),alpha(s(1:(S_steps/L):end)),s(1:end));
gamma = (((s-1)*Ds).*(L-(s-1)*Ds)/L^2)';
%f = @(time) (3-L)/6*(sin((time-1)*Dt)+5);
f = @(time) 0;
lambda1 = 0;
lambda2 = 1 - lambda1;
phi = 1; % потенциально функция фи в 1 функционале
p = 1; % коэффициент в 1 функционале
%eps = 1/100;
eps = 0;
umax = 10000;
umin = 0;
%% Initialization
CorrectX(s(:), t(:)) = ((sin((t-1)*Dt) + 5)'*(6*(s-1)*Ds/L^2+1))';
CorrectU(s(:), t(:)) = ((sin((t-1)*Dt) + 5)'*((alpha(s)-1).*(6*(s-1)*Ds/L^2 + 1) - 6/L^2) - cos((t-1)*Dt)'*(6*(s-1)*Ds/L^2 + 1))';
%% Заполняем краевые условия
CorrectU(s(:), t(:)) = zeros(size(s,2), size(t,2));
temp = [2 2 2 2 620 2 8 2 1603; 9340 207 524 1431 982 1794 108 767 5005; 7691 5962 2151 1416 2211 4492 1275 4420 5740; 1959 2149 1206 1088 675 3841 1779 6569 3534; 887 1395 1986 5585 576 1914 2442 8441 3473; 474 413 1401 2610 269 917 1465 5248 1829; 268 134 195 142 104 332 679 2535 748; 2 52 30 2 38 2 405 200 284];
CorrectU(s(1:(S_steps/L):end), t(1:(T_steps/T):end)) = temp;
for j=t
    CorrectU(s, j) = interp1(s(1:(S_steps/L):end),CorrectU(s(1:(S_steps/L):end), j),s(1:end));
end
for i = s
    CorrectU(i, t) = interp1(t(1:(T_steps/T):end),CorrectU(i, t(1:(T_steps/T):end)),t(1:end));
end
CorrectU(s(:), t(:)) = zeros(size(s,2), size(t,2));
CorrectX(s(1:(S_steps/L):end), 1) = [280026 118907 33972 10172 2456 993 483 3]';
CorrectX(s, 1) = interp1(s(1:(S_steps/L):end),CorrectX(s(1:(S_steps/L):end), 1),s(1:end));

%% Iteration Scheme
xu = Boundary(CorrectX(s,1), CorrectU);

%% Создание х с волной для второго функционала
need_x(s) = xu(s, end);
%need_x(s(1:(S_steps/L):end)) = [813983 150629 52499 20628 9938 3524 986 374]';
%need_x(s) = interp1(s(1:(S_steps/L):end), need_x(s(1:(S_steps/L):end)), s(1:end));

%% Подсчет функционала для аналитического решения
J2uCorrect = trapz(0:Ds:L, (xu(s, end) - need_x(s)').^2)
J1uCorrect = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectU(s, t)))
%% всякие графики с разными управлениями
% неплохо потенциально сделать с базисными
%{
for i=-5:0.1:5
    CorrectUh = i*ones(size(s,2), size(t,2));
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, 'bx', 'markers', 10);
end
% u(s,t) = s + i
for i=-5:0.1:5
    for time=t
        CorrectUh(:,time) = (s-1)*Ds + i;
    end
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, 'bo', 'markers', 10);
end
% u(s,t) = t + i
for i=-5:0.1:5
    for x=s
        CorrectUh(x,:) = (t-1)*Dt + i;
    end
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, 'ro', 'markers', 10);
end   
%u(s,t) = t*s+i;
for i=-5:0.1:5
    for x=s
        for time=t
            CorrectUh(x,time) = (time-1)*Dt * (x-1)*Ds + i;
        end
    end
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, 'yo', 'markers', 10);
end   
%u(s,t) = t+s+i;
for i=-5:0.1:5
    for x=s
        for time=t
            CorrectUh(x,time) = (time-1)*Dt + (x-1)*Ds + i;
        end
    end
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, 'ko', 'markers', 10);
end   
%u(s,t) = t-s+i;
for i=-5:0.1:5
    for x=s
        for time=t
            CorrectUh(x,time) = (time-1)*Dt - (x-1)*Ds + i;
        end
    end
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, 'mo', 'markers', 10);
end   
%u(s,t) = s-t+i;
for i=-5:0.1:5
    for x=s
        for time=t
            CorrectUh(x,time) = - (time-1)*Dt + (x-1)*Ds + i;
        end
    end
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, 'go', 'markers', 10);
end   
for i=-5:0.1:5
    CorrectUh = i*ones(size(s,2), size(t,2));
    CorrectUh(1:end, end/2:end) = (i+1);
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, 'yx', 'markers', 10);
end
for i=-5:0.1:5
    CorrectUh = i*ones(size(s,2), size(t,2));
    CorrectUh(1:end, end/2:end) = -i;
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, 'mx', 'markers', 10);
end
for i=-5:0.1:5
    CorrectUh = (i+1)*ones(size(s,2), size(t,2));
    CorrectUh(1:end, end/2:end) = i;
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, 'gx', 'markers', 10);
end
for i=-5:0.1:5
    CorrectUh = i*ones(size(s,2), size(t,2));
    CorrectUh(1:end, end/3:2*end/3) = (i+1);
    CorrectUh(1:end, 2*end/3:end) = (i+2);
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, 'rx', 'markers', 10);
end
for i=1:100
    CorrectUh = -10*rand(size(s,2),size(t,2))+5;
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, 'o');
end

for i=1:100
    CorrectUh = CorrectU+10*rand(size(s,2),size(t,2));
    xuh(s, 1) = CorrectX(s, 1);
    for time=t(1:end - 1)
    
        xuh(1, time) = 1/2*trapz(xuh(2:end, time))*Ds - 5/4*L*sin((time-1)*Dt) - 5/4*L + 2*(sin((time-1)*Dt)+1);    
    
        xuh(s(2:end), time + 1) = xuh(s(2:end), time) + Dt*((alpha-1)*xuh(s(2:end)-1, time) - CorrectUh(s(2:end)-1, time)) - (Dt/Ds)*(xuh(s(2:end), time)-xuh(s(2:end)-1, time));
        
    end
    J2uh = trapz(0:Ds:L, (xuh(s, end) - need_x(s)').^2);
    J1uh = p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectUh(s, t)));
    plot(J2uh, -J1uh, '*');
end
%}
%{
figure
hold on;
plot(CorrectX(end,:),'-')
plot(xu(end, :),'--')
legend({'Аналитическое Решение','Численное Решение'});
set(gca, 'XtickLabel', {0:0.5:4});
xlabel('t')
ylabel('x(L, t)')
figure
hold on;
plot(CorrectX(:,t(end)),'-')
plot(xu(:, t(end)),'--')
legend({'Аналитическое Решение','Численное Решение'});
set(gca, 'XtickLabel', {0:0.5:3});
xlabel('s')
ylabel('x(s, T)')
figure
hold on;
plot(CorrectX(:,1*T_steps/T+1),'-')
plot(xu(:, 1*T_steps/T+1),'--')
legend({'Аналитическое Решение','Численное Решение'});
set(gca, 'XtickLabel', {0:0.5:4});
xlabel('s')
ylabel('x(s, 1)')
figure
hold on;
plot(CorrectX(:,2*T_steps/T+1),'-')
plot(xu(:, 2*T_steps/T+1),'--')
legend({'Аналитическое Решение','Численное Решение'});
set(gca, 'XtickLabel', {0:0.5:4});
xlabel('s')
ylabel('x(s, 2)')
figure
hold on;
plot(CorrectX(1,:),'-')
plot(xu(1, :),'--')
legend({'Аналитическое Решение','Численное Решение'});
set(gca, 'XtickLabel', {0:0.5:4});
xlabel('t')
ylabel('x(0, t)')
hold off;
%}
%% Цикл решения задачи экстраградиентным методом

% создание матрицы новых управлений

figure(1)
hold on;
lambda1 = 0;
while (lambda1 <= 1)
    k = 1;
    J2u = 1;
    uk = 5000*ones(size(s,2), size(t,2));
    uk_ = zeros(size(s,2), size(t,2));
    uprev = -1*ones(size(s,2), size(t,2));
    prev = 1;
while (J2u > 10e-3) && (k < 1000) && (prev > 10e-3)
    % Прямая задача
    xu = Boundary(CorrectX(s,1), uk);

    % Обратная краевая задача
    psiStartCondition = -2*(xu(s, end) - need_x(s)');
    psi = reverseBoundary(psiStartCondition);
    
    %% Подсчет прогнозного шага
    % Подсчет производных функционалов
    J1_ = -p*phi*L*T;
    J2_ = psi;
    
    beta = (500+k)*sqrt(Ds^2+Dt^2)/norm(J2_);
    if beta < 0.1
        beta = 0.1;
    end
    
    % Проекция
    uk_ = projection(uk - beta * (lambda1 * J1_ + lambda2 * J2_));
    
    %% Подсчет основного шага
    % Прямая задача
    xu = Boundary(CorrectX(s,1), uk_);

    % Обратная задача
    psiStartCondition = -2*(xu(s, end) - need_x(s)');
    psi = reverseBoundary(psiStartCondition);
    
    % Подсчет производных функционалов
    J1_ = -p*phi*L*T;
    J2_ = psi;
    
    beta = (500+k)*sqrt(Ds^2+Dt^2)/norm(J2_);
    if beta < 0.1
        beta = 0.1;
    end
    
    % Проекция
    uk = projection(uk_ - beta * (lambda1 * J1_ + lambda2 * J2_));
    
    %Вычисление решения
    xu = Boundary(CorrectX(s, 1), uk);
    
    k
    % Подсчет функционала для численного решения
    J2u = trapz((xu(s, end) - need_x(s)').^2)*Ds
    J1u = -p*trapz(trapz(phi*uk(s, t))*Ds)*Dt
    k = k+1;
    prev = sum(sum((uk-uprev).^2))
    uprev = uk;
end
    figure(1)
    if (lambda1 == 0)
        plot(J2u, J1u, 'ro', 'markers', 5);
    elseif (lambda1 == 1)
        plot(J2u, J1u, 'yo', 'markers', 5);
    else
        plot(J2u, J1u, 'o', 'markers', 5, 'Color', [0, lambda1, 0] );
    end
    figure
    hold on;
    plot(CorrectX(end,:),'-')
    plot(xu(end, :),'--')
    legend({'Аналитическое Решение','Численное Решение'});
    set(gca, 'XtickLabel', {0:0.5:4});
    xlabel('t')
    ylabel('x(L, t)')
    figure
    hold on;
    plot(CorrectX(:,t(end)),'-')
    plot(xu(:, t(end)),'--')
    legend({'Аналитическое Решение','Численное Решение'});
    set(gca, 'XtickLabel', {0:0.5:3});
    xlabel('s')
    ylabel('x(s, T)')
    figure
    hold on;
    plot(CorrectX(:,1*T_steps/T+1),'-')
    plot(xu(:, 1*T_steps/T+1),'--')
    legend({'Аналитическое Решение','Численное Решение'});
    set(gca, 'XtickLabel', {0:0.5:4});
    xlabel('s')
    ylabel('x(s, 1)')
    figure
    hold on;
    plot(CorrectX(:,2*T_steps/T+1),'-')
    plot(xu(:, 2*T_steps/T+1),'--')
    legend({'Аналитическое Решение','Численное Решение'});
    set(gca, 'XtickLabel', {0:0.5:4});
    xlabel('s')
    ylabel('x(s, 2)')
    figure
    hold on;
    plot(CorrectX(1,:),'-')
    plot(xu(1, :),'--')
    legend({'Аналитическое Решение','Численное Решение'});
    set(gca, 'XtickLabel', {0:0.5:4});
    xlabel('t')
    ylabel('x(0, t)')
    hold off;
    figure
    hold on;
    plot(CorrectU(end,:),'-')
    plot(uk(end, :),'--')
    legend({'Аналитическое Решение','Численное Решение'});
    set(gca, 'XtickLabel', {0:0.5:4});
    xlabel('t')
    ylabel('u(L, t)')
    figure
    hold on;
    plot(CorrectU(:,t(end)),'-')
    plot(uk(:, t(end)),'--')
    legend({'Аналитическое Решение','Численное Решение'});
    set(gca, 'XtickLabel', {0:0.5:3});
    xlabel('s')
    ylabel('u(s, T)')
    figure
    hold on;
    plot(CorrectU(:,1*T_steps/T+1),'-')
    plot(uk(:, 1*T_steps/T+1),'--')
    legend({'Аналитическое Решение','Численное Решение'});
    set(gca, 'XtickLabel', {0:0.5:4});
    xlabel('s')
    ylabel('u(s, 1)')
    figure
    hold on;
    plot(CorrectU(:,2*T_steps/T+1),'-')
    plot(uk(:, 2*T_steps/T+1),'--')
    legend({'Аналитическое Решение','Численное Решение'});
    set(gca, 'XtickLabel', {0:0.5:4});
    xlabel('s')
    ylabel('u(s, 2)')
    figure
    hold on;
    plot(CorrectU(1,:),'-')
    plot(uk(1, :),'--')
    legend({'Аналитическое Решение','Численное Решение'});
    set(gca, 'XtickLabel', {0:0.5:4});
    xlabel('t')
    ylabel('u(0, t)')
    hold off;
    lambda1 = lambda1+2;
    lambda2 = 1 - lambda1;
end
xlabel('J2u');
ylabel('-J1u');
hold off;