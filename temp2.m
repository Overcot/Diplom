clc; clear;

%% Params of Setki
global Ds Dt s t alpha f gamma precision;
precision = 1;
folder_to_save = ['tahvonen&ours(1 & ',num2str(precision), ')'];

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
format long
beta = 11.269799;
while (abs(CorrectX(1,1) - beta*trapz(gamma(1:end).*CorrectX(1:end, 1))*Ds) >= 10e-9)
    beta = beta + 10e-12;
    abs(CorrectX(1,1) - beta*trapz(gamma(1:end).*CorrectX(1:end, 1))*Ds)
end
beta
CorrectX(1,1) - beta*trapz(gamma(1:end).*CorrectX(1:end, 1)*Ds)