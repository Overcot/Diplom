clc; clear;

%% Params of Setki
global Ds Dt s t alpha f gamma umin umax;
folder_to_save = 'data_when_u=0(2)';

L = 7;
T = 8;
Ds = 0.01;
Dt = 0.01;
S_steps= L/Ds;
T_steps = T/Dt;
s=[1:S_steps+1];
t=[1:T_steps+1];
alpha(s(1:(S_steps/L):end)) = [1 0.95 0.85 0.8 0.7 0.5 0.3 0];
alpha(s) = interp1(s(1:(S_steps/L):end),alpha(s(1:(S_steps/L):end)),s(1:end));
gamma = (((s-1)*Ds).*(L-(s-1)*Ds)/L^2)';

f = @(time) 0;

umax = 10000;
umin = 0;
%% Initialization

CorrectU(s(:), t(:)) = zeros(size(s,2), size(t,2));
u_input = [2 2 2 2 620 2 8 2 1603; 9340 207 524 1431 982 1794 108 767 5005; 7691 5962 2151 1416 2211 4492 1275 4420 5740; 1959 2149 1206 1088 675 3841 1779 6569 3534; 887 1395 1986 5585 576 1914 2442 8441 3473; 474 413 1401 2610 269 917 1465 5248 1829; 268 134 195 142 104 332 679 2535 748; 2 52 30 2 38 2 405 200 284];
CorrectU(s(1:(S_steps/L):end), t(1:(T_steps/T):end)) = u_input;
for j=t
    CorrectU(s, j) = interp1(s(1:(S_steps/L):end),CorrectU(s(1:(S_steps/L):end), j),s(1:end));
end
for i = s
    CorrectU(i, t) = interp1(t(1:(T_steps/T):end),CorrectU(i, t(1:(T_steps/T):end)),t(1:end));
end

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
%% Iteration Scheme
xu = Boundary(CorrectX(s,1), CorrectU);
    mkdir(folder_to_save); 
    plotGraph(CorrectX(1,:), xu(1, :), {0:T}, 't', 'x(s=0, t)', folder_to_save);
    plotGraph(CorrectX(2*S_steps/L+1,:), xu(2*S_steps/L+1, :), {0:T}, 't', 'x(s=2, t)', folder_to_save);
    plotGraph(CorrectX(5*S_steps/L+1,:), xu(5*S_steps/L+1, :), {0:T}, 't', 'x(s=5, t)', folder_to_save);
    plotGraph(CorrectX(end,:), xu(end, :), {0:T}, 't', 'x(s=7, t)', folder_to_save);
    plotGraph(CorrectX(:,1), xu(:,1), {0:L}, 's', 'x(s, t=0)', folder_to_save);
    plotGraph(CorrectX(:,2*T_steps/T+1), xu(:, 2*T_steps/T+1), {0:L}, 's', 'x(s, t=2)', folder_to_save);
    plotGraph(CorrectX(:,5*T_steps/T+1), xu(:, 5*T_steps/T+1), {0:L}, 's', 'x(s, t=5)', folder_to_save);
    plotGraph(CorrectX(:,end), xu(:,end), {0:L}, 's', 'x(s, t=8)', folder_to_save);