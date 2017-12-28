clc; clear;

%% Params of Setki
global Ds Dt s t alpha f gamma umin umax precision;
precision = 0.01;


folder_to_save = ['0_control(',num2str(precision), ')A&B'];

L = 7;
T = 8;
Ds = precision;
Dt = precision;
S_steps= L/Ds;
T_steps = T/Dt;
s=1:S_steps+1;
t=1:T_steps+1;

alpha(s(1:(S_steps/L):end)) = [1 0.95 0.85 0.8 0.7 0.5 0.3 0];
alpha(s) = interp1(s(1:(S_steps/L):end),alpha(s(1:(S_steps/L):end)),s(1:end));

f = @(time) 0;


%% Init
x_input = [280026 558880 682092 670578 557756 537638 691972 753642 813983; 
    118907 55968 111702 136328 134027 111239 107457 138303 150629; 
    33972 40021 21307 42458 51348 50734 41526 41080 52499; 
    10172 14084 18868 10804 23676 28260 26195 23240 20628; 
    2456 5064 7454 11312 6161 14862 15329 15621 9938; 
    993 879 2141 3203 2933 3475 7971 7856 3524; 
    483 254 229 288 127 1590 1429 3760 986; 
    3 83 48 3 61 3 649 320 374];
CorrectX = input_data(x_input, S_steps, T_steps, L, T);

%A = (trapz(gamma(1:end).*CorrectX(1:end, 1))*Ds - CorrectX(1,1))/(1-trapz(gamma(1:end))*Ds)
%CorrectX(s,1) = CorrectX(s,1) + A;
%% Iteration Scheme
CorrectU(s(:), t(:)) = zeros(size(s,2), size(t,2));

if precision == 1
    beta_for_gamma = 11.269799576773186;
elseif precision == 0.1
    beta_for_gamma = 9.207024636678527;
elseif precision == 0.01
    beta_for_gamma = 9.190203289586552;
else
    %precision == 0.001
    beta_for_gamma = 9.190035386505700;
end
gamma = beta_for_gamma*(((s-1)*Ds).*(L-(s-1)*Ds)/L^2)';
xuB(s(:), t(:)) = zeros(size(s,2), size(t,2));
xuB = Boundary(CorrectX(s,1), CorrectU);


beta_for_gamma = 1;
gamma = beta_for_gamma*(((s-1)*Ds).*(L-(s-1)*Ds)/L^2)';
A = (trapz(gamma(1:end).*CorrectX(1:end, 1))*Ds - CorrectX(1,1))/(1-trapz(gamma(1:end))*Ds)
CorrectX(s,1) = CorrectX(s,1) + A;
xuA(s(:), t(:)) = zeros(size(s,2), size(t,2));
xuA = Boundary(CorrectX(s,1), CorrectU);


mkdir(folder_to_save);
    plotGraph2(xuB(1, :), xuA(1, :), {0:T}, 't', 'x(s=0, t)', folder_to_save, 'beta', 'A');
    plotGraph2(xuB(2*S_steps/L+1, :), xuA(2*S_steps/L+1, :), {0:T}, 't', 'x(s=2, t)', folder_to_save, 'beta', 'A');
    plotGraph2(xuB(5*S_steps/L+1, :), xuA(5*S_steps/L+1, :), {0:T}, 't', 'x(s=5, t)', folder_to_save, 'beta', 'A');
    plotGraph2(xuB(end, :), xuA(end, :), {0:T}, 't', 'x(s=7, t)', folder_to_save, 'beta', 'A');
    plotGraph2(xuB(:,1), xuA(:,1), {0:L}, 's', 'x(s, t=0)', folder_to_save, 'beta', 'A');
    plotGraph2(xuB(:,2*T_steps/T+1), xuA(:, 2*T_steps/T+1), {0:L}, 's', 'x(s, t=2)', folder_to_save, 'beta', 'A');
    plotGraph2(xuB(:,5*T_steps/T+1), xuA(:, 5*T_steps/T+1), {0:L}, 's', 'x(s, t=5)', folder_to_save, 'beta', 'A');
    plotGraph2(xuB(:,end), xuA(:,end), {0:L}, 's', 'x(s, t=8)', folder_to_save, 'beta', 'A');
