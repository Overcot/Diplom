clc; clear;

global Ds Dt s t alpha gamma precision;
precision = 0.01;
folder_to_save = ['model(',num2str(precision), ')'];
L = 10; % because in data they start from 1 not 0
T = 11;
Ds = precision;
Dt = precision;
S_steps= L/Ds;
T_steps = T/Dt;
s=[1:S_steps+1];
t=[1:T_steps+1];

alpha(s(1:(S_steps/L):end)) = [0.319 0.409 0.788 0.818 0.818 0.818 0.818 0.818 0.818 0.818 0.818 ];
alpha(s) = interp1(s(1:(S_steps/L):end),alpha(s(1:(S_steps/L):end)),s(1:end));
gamma(s(1:(S_steps/L):end)) = [0.002728 0.098784 0.7585 2.643219 5.34322 8.227 9.891 11.075 12.397 13.494 15.016];
gamma(s) = interp1(s(1:(S_steps/L):end),gamma(s(1:(S_steps/L):end)),s(1:end));
%%

CorrectU(s(:), t(:)) = zeros(size(s,2), size(t,2));

CorrectX(s(:), t(:)) = zeros(size(s,2), size(t,2));
x_input = [280026 558880 682092 670578 557756 537638 691972 753642 813983; 118907 55968 111702 136328 134027 111239 107457 138303 150629; 33972 40021 21307 42458 51348 50734 41526 41080 52499; 10172 14084 18868 10804 23676 28260 26195 23240 20628; 2456 5064 7454 11312 6161 14862 15329 15621 9938; 993 879 2141 3203 2933 3475 7971 7856 3524; 483 254 229 288 127 1590 1429 3760 986; 3 83 48 3 61 3 649 320 374];
%CorrectX(s(301:(S_steps/L):end), t(1:(T_steps/T):end)) = x_input;
for j=t
    CorrectX(s, j) = interp1(s(1:(S_steps/L):end),CorrectX(s(1:(S_steps/L):end), j),s(1:end));
end
for i = s
    CorrectX(i, t) = interp1(t(1:(T_steps/T):end),CorrectX(i, t(1:(T_steps/T):end)),t(1:end));
end

%% Our model
xu2(s(:), t(:)) = zeros(size(s,2), size(t,2));
xu2(s(1:(S_steps/L):end-300), 1) = x_input(:,1);
xu2(s, 1) = interp1(s(1:(S_steps/L):end),xu2(s(1:(S_steps/L):end), 1),s(1:end));
u(s(:), t(:)) = zeros(size(s,2), size(t,2));


xu2 = Boundary(xu2(s, 1), CorrectU);

%% Graphs
mkdir(folder_to_save);
2*S_steps/L+1
    plotGraph2(xu2(1, :), CorrectX(1, :), {0:T}, 't', 'x(s=0, t)', folder_to_save);
    plotGraph2(xu2(2*S_steps/L+1, :), CorrectX(2*S_steps/L+1, :), {0:T}, 't', 'x(s=2, t)', folder_to_save);
    plotGraph2(xu2(5*S_steps/L+1, :), CorrectX(5*S_steps/L+1, :), {0:T}, 't', 'x(s=5, t)', folder_to_save);
    plotGraph2(xu2(end, :), CorrectX(end, :), {0:T}, 't', 'x(s=7, t)', folder_to_save);
    plotGraph2(xu2(:,1), CorrectX(:,1), {0:L}, 's', 'x(s, t=0)', folder_to_save);
    plotGraph2(xu2(:,2*T_steps/T+1), CorrectX(:, 2*T_steps/T+1), {0:L}, 's', 'x(s, t=2)', folder_to_save);
    plotGraph2(xu2(:,5*T_steps/T+1), CorrectX(:, 5*T_steps/T+1), {0:L}, 's', 'x(s, t=5)', folder_to_save);
    plotGraph2(xu2(:,end), CorrectX(:,end), {0:L}, 's', 'x(s, t=8)', folder_to_save);