clc; clear;

%% Params of Setki
global Ds Dt s t alpha f gamma;
folder_to_save = 'tahvonen&ours(1 & 0.1)';

L = 7;
T = 8;
Ds = 1;
Dt = 1;
S_steps_tahn= L/Ds;
T_steps_tahn = T/Dt;
s=[1:S_steps_tahn+1];
t=[1:T_steps_tahn+1];
alpha(s(1:(S_steps_tahn/L):end)) = [1 0.95 0.85 0.8 0.7 0.5 0.3 0];
alpha(s) = interp1(s(1:(S_steps_tahn/L):end),alpha(s(1:(S_steps_tahn/L):end)),s(1:end));
gamma = (((s-1)*Ds).*(L-(s-1)*Ds)/L^2)';

f = @(time) 0;


%% Initialization

CorrectU(s(:), t(:)) = zeros(size(s,2), size(t,2));

CorrectX(s(:), t(:)) = zeros(size(s,2), size(t,2));
x_input = [280026 558880 682092 670578 557756 537638 691972 753642 813983; 118907 55968 111702 136328 134027 111239 107457 138303 150629; 33972 40021 21307 42458 51348 50734 41526 41080 52499; 10172 14084 18868 10804 23676 28260 26195 23240 20628; 2456 5064 7454 11312 6161 14862 15329 15621 9938; 993 879 2141 3203 2933 3475 7971 7856 3524; 483 254 229 288 127 1590 1429 3760 986; 3 83 48 3 61 3 649 320 374];
CorrectX(s(1:(S_steps_tahn/L):end), t(1:(T_steps_tahn/T):end)) = x_input;
for j=t
    CorrectX(s, j) = interp1(s(1:(S_steps_tahn/L):end),CorrectX(s(1:(S_steps_tahn/L):end), j),s(1:end));
end
for i = s
    CorrectX(i, t) = interp1(t(1:(T_steps_tahn/T):end),CorrectX(i, t(1:(T_steps_tahn/T):end)),t(1:end));
end
%% Tahnoven
xu(s(:), t(:)) = zeros(size(s,2), size(t,2));
u(s(:), t(:)) = zeros(size(s,2), size(t,2));
xu(s, 1) = CorrectX(s, 1);

for time=t(1:end - 1)
    xu(1, time) = trapz(gamma(1:end).*xu(1:end, time))*Ds + f(time);
    xu(s(2:end), time + 1) = alpha(s(2:end)-1)'.*xu(s(2:end)-1, time) - u(s(2:end)-1, time);
end
xu(1, end) = trapz(gamma(1:end).*xu(1:end, t(end)))*Ds + f(t(end));

%% Our model

Ds = 0.1;
Dt = 0.1;
S_steps= L/Ds;
T_steps = T/Dt;
s=[1:S_steps+1];
t=[1:T_steps+1];
alpha(s(1:(S_steps/L):end)) = [1 0.95 0.85 0.8 0.7 0.5 0.3 0];
alpha(s) = interp1(s(1:(S_steps/L):end),alpha(s(1:(S_steps/L):end)),s(1:end));
gamma = (((s-1)*Ds).*(L-(s-1)*Ds)/L^2)';

%%

CorrectU(s(:), t(:)) = zeros(size(s,2), size(t,2));

CorrectX(s(:), t(:)) = zeros(size(s,2), size(t,2));
x_input = [280026 558880 682092 670578 557756 537638 691972 753642 813983; 118907 55968 111702 136328 134027 111239 107457 138303 150629; 33972 40021 21307 42458 51348 50734 41526 41080 52499; 10172 14084 18868 10804 23676 28260 26195 23240 20628; 2456 5064 7454 11312 6161 14862 15329 15621 9938; 993 879 2141 3203 2933 3475 7971 7856 3524; 483 254 229 288 127 1590 1429 3760 986; 3 83 48 3 61 3 649 320 374];
CorrectX(s(1:(S_steps/L):end), t(1:(T_steps/T):end)) = x_input;
for j=t
    CorrectX(s, j) = interp1(s(1:(S_steps/L):end),CorrectX(s(1:(S_steps/L):end), j),s(1:end));
end
for i = s
    CorrectX(i, t) = interp1(t(1:(T_steps/T):end),CorrectX(i, t(1:(T_steps/T):end)),t(1:end));
end

%% Our model
xu2(s(:), t(:)) = zeros(size(s,2), size(t,2));
u(s(:), t(:)) = zeros(size(s,2), size(t,2));
xu2(s, 1) = CorrectX(s, 1);

xu2 = Boundary(CorrectX(s,1), CorrectU);

%% Graphs
    mkdir(folder_to_save);
    

temp = xu;
xu = zeros(size(s,2), size(t,2));
xu(s(1:(S_steps/L):end), t(1:(T_steps/T):end)) = temp;
for j=t
    xu(s, j) = interp1(s(1:(S_steps/L):end),xu(s(1:(S_steps/L):end), j),s(1:end));
end
for i = s
    xu(i, t) = interp1(t(1:(T_steps/T):end),xu(i, t(1:(T_steps/T):end)),t(1:end));
end

    plotGraph2(xu2(1, :), xu(1, :), {0:T}, 't', 'x(s=0, t)', folder_to_save);
    plotGraph2(xu2(2*S_steps/L+1, :), xu(2*S_steps/L+1, :), {0:T}, 't', 'x(s=2, t)', folder_to_save);
    plotGraph2(xu2(5*S_steps/L+1, :), xu(5*S_steps/L+1, :), {0:T}, 't', 'x(s=5, t)', folder_to_save);
    plotGraph2(xu2(end, :), xu(end, :), {0:T}, 't', 'x(s=7, t)', folder_to_save);
    plotGraph2(xu2(:,1), xu(:,1), {0:L}, 's', 'x(s, t=0)', folder_to_save);
    plotGraph2(xu2(:,2*T_steps/T+1), xu(:, 2*T_steps/T+1), {0:L}, 's', 'x(s, t=2)', folder_to_save);
    plotGraph2(xu2(:,5*T_steps/T+1), xu(:, 5*T_steps/T+1), {0:L}, 's', 'x(s, t=5)', folder_to_save);
    plotGraph2(xu2(:,end), xu(:,end), {0:L}, 's', 'x(s, t=8)', folder_to_save);