clear; clc;
global ssb a b
precision = 0.1;
L = 5;
T = 5;
gammaValue = 1/L;
muValue = 1/2;
x0 = 1;
%anna model
a = 8.244595391;
b = 8.58578E-07;
allee = 0.2;
%% Numerical Part
[x, L, T, Ds, Dt, S_steps, T_steps, s, t, gamma, mu, p] = startInit(precision, L, T, gammaValue, muValue, x0);

for time=t(1:end - 1)
    
    x(1, time) = recruitmentFunction('Tahvonen', time, x, [Ds, Dt, gamma, p(time)]);
    %x(1, time) = recruitmentFunction('Anna', time, x, [Ds, Dt, gamma, p(time), a, b, allee]);
    %%
    %{
    %{
    for class=s(1:end-1)
        x(class, time+1) = x(class, time) - Dt*(mu(class)*x(class, time) + (x(class+1, time) - x(class, time))/Ds);
    end
    %}
    %{
    for class = s(1:end-1)
        x(class+1, time+1) = x(class+1, time) - Dt*((x(class+1, time)-x(class, time))/(Ds) + mu(class+1)*x(class+1, time));
    end
    %}
    %{
    for class = s(2:end-1)
        x(class, time+1) = x(class, time) - Dt*((x(class+1, time)-x(class-1, time))/(2*Ds) + mu(class)*x(class, time));
    end
    %}
    %}
    %% Potapov
    for class = s(1:end-1)
        %(x(class+1, time+1) - x(class, time))/Dt = -mu(class) * (x(class+1, time+1)+x(class, time))/2;
        %x(class+1,time+1)/Dt + mu(class)*x(class+1, time+1)/2 = x(class, time)/Dt - mu(class)*x(class, time)/2;
        x(class+1, time+1) = (x(class, time)/Dt - mu*x(class, time)/2)/(1/Dt + mu/2);
    end
end
x(1,t(end)) = (p(t(end))) + trapz(gamma*x(2:end, t(end))*Ds)/(1-gamma*Ds);

     
%% Analytical Check
t_new = 1:S_steps+1;

solution = zeros(size(s,2),size(s,2));

for i = 1:T/L
    for time=(i-1)*S_steps+1:i*S_steps+1
        %(i-1)*S_steps*Ds - (time-1)*Dt %= KL-t
        %i*S_steps+1 - (time-1)
        %gamma*exp(mu*((i-1)*S_steps*Ds - (time-1)*Dt) *
        %trapz(x(1 : i*S_steps*Ds - (time-1)*Dt), (i-1)*S_steps+1))*Ds
        ksi(time) = p(time) + gamma*exp(mu*((i-1)*S_steps*Ds - (time-1)*Dt))*trapz(x(1 : i*S_steps+1 - (time-1), (i-1)*S_steps+1)*Ds);
    end
    a(1) = ksi(1);
    for time=(i-1)*S_steps+1 + 1:i*S_steps+1 % add 1 more because we need 2 points at start
        k = ((i-1)*S_steps:time-1)*Dt; % 0 0.1 0.2 0.3
        l = ((i-1)*S_steps+1:time);
        %exp((mu-gamma)*k).*ksi(l)*Dt;
        a(time) = ksi(time)+gamma*exp((gamma-mu)*(time-1)*Dt)*trapz(exp((mu-gamma)*k).*ksi(l)*Dt);
    end
    for time=(i-1)*S_steps+1:i*S_steps+1
        for class=s(1:end)
            if class < time - (i-1)*S_steps+1
                solution(class, time) = a(time-class+1)*exp(-mu*(class-1)*Ds);
            else
                %exp((i-1)*S_steps*Ds - (time-1)*Dt+class);
                %(i-1)*S_steps+1-time+class;
                %x((i-1)*S_steps+1-time+class, (i-1)*S_steps+1);
                solution(class, time) =x((i-1)*S_steps+1-time+class, (i-1)*S_steps+1) * exp(mu*((i-1)*S_steps*Ds - (time-1)*Dt));
            end
        end
    end
end


sum_of_differences = 0;
for time = t(1:end)
    for class=s(1:end)
        sum_of_differences = sum_of_differences + abs(x(class, time) - solution(class, time));
    end
end
sum_of_differences
otnositelnaya = trapz(trapz((x-solution).^2*Ds)*Dt)/trapz(trapz((solution.^2)*Ds)*Dt)*100

%{
figure
hold on;
plot(diag(x))
plot(diag(solution), '-')
legend('numerical','analytical')
title('numerical and analytical solution x(t,t)')


figure
hold on;
plot(diag(x, -1/precision))
plot(diag(solution, -1/precision))
legend('numerical','analytical')
title('numerical and analytical solution x(t+1,t)')


figure
hold on;
plot(diag(x,2/precision))
plot(diag(solution,2/precision))
legend('numerical','analytical')
title('numerical and analytical solution x(t-2,t)')
%}
figure
hold on;
plot(x(1,t))
plot(solution(1,t),'o','MarkerSize',2)

xlim([0 T_steps+1])

set(gca, 'XtickLabel', {0:1:T+1});
xlabel('t')

legend('numerical','analytical')
title('x(0,t) & a(t)')


figure
hold on;
plot(x(S_steps/2+1,t))
plot(solution(S_steps/2+1,t),'o','MarkerSize',2)

xlim([0 T_steps+1])

set(gca, 'XtickLabel', {0:1:T+1});
xlabel('t')

legend('numerical','analytical')
title('x(L/2,t) & solution(L/2,t)')

figure
hold on;
plot(x(S_steps+1,t))
plot(solution(S_steps+1,t),'o','MarkerSize',2)

xlim([0 T_steps+1])

set(gca, 'XtickLabel', {0:1:T+1});
xlabel('t')

legend('numerical','analytical')
title('x(L,t) & solution(L,t)')



figure
hold on;
plot(x(s,T_steps/2+1))
plot(solution(s,T_steps/2+1),'o','MarkerSize',2)

xlim([0 S_steps+1])

set(gca, 'XtickLabel', {0:L+1});
xlabel('s')

legend('numerical','analytical')
title('x(s, T/2) & solution(s, T/2)')

figure
hold on;
plot(x(s,T_steps+1))
plot(solution(s,T_steps+1),'o','MarkerSize',2)

xlim([0 S_steps+1])

set(gca, 'XtickLabel', {0:L+1});
xlabel('s')

legend('numerical','analytical')
title('x(s,T) & solution(s,T)')

%}



%{
string_numerical = ['numerical = ',num2str(x(1,1+5/precision))];
text(1+5/precision+50,x(1,1+5/precision),string_numerical,'FontSize', 12);
string_analytical = ['analytical = ',num2str(solution(1,1+5/precision))];
text(1+5/precision+50,solution(1,1+5/precision),string_analytical,'FontSize', 12);
line([1+5/precision 1+5/precision], [0 solution(1, 1+5/precision)], 'LineStyle','--');

string_numerical = ['numerical = ',num2str(x(1,1+10/precision))];
text(1+10/precision+20,x(1,1+10/precision),string_numerical,'FontSize', 12);
string_analytical = ['analytical = ',num2str(solution(1,1+10/precision))];
text(1+10/precision+20,solution(1,1+10/precision),string_analytical,'FontSize', 12);
line([1+10/precision 1+10/precision], [0 solution(1, 1+10/precision)], 'LineStyle','--');
%}



