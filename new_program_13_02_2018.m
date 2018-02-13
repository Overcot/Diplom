clc; clear;

%% Params of Setki
global Ds Dt s t alpha f gamma umin umax;
precision = 0.01;

folder_to_save = ['new_prog(',num2str(precision), ')'];

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

lambda1 = 0;
lambda2 = 1 - lambda1;
phi = 1;
p = 1;

eps = 0;
umax = 10000;
umin = 0;

%% Init
u_input = [2 2 2 2 620 2 8 2 1603; 
    9340 207 524 1431 982 1794 108 767 5005; 
    7691 5962 2151 1416 2211 4492 1275 4420 5740; 
    1959 2149 1206 1088 675 3841 1779 6569 3534; 
    887 1395 1986 5585 576 1914 2442 8441 3473; 
    474 413 1401 2610 269 917 1465 5248 1829; 
    268 134 195 142 104 332 679 2535 748; 
    2 52 30 2 38 2 405 200 284];

CorrectU = input_data(u_input, S_steps, T_steps, L, T);

x_input = [280026 558880 682092 670578 557756 537638 691972 753642 813983; 
    118907 55968 111702 136328 134027 111239 107457 138303 150629; 
    33972 40021 21307 42458 51348 50734 41526 41080 52499; 
    10172 14084 18868 10804 23676 28260 26195 23240 20628; 
    2456 5064 7454 11312 6161 14862 15329 15621 9938; 
    993 879 2141 3203 2933 3475 7971 7856 3524; 
    483 254 229 288 127 1590 1429 3760 986; 
    3 83 48 3 61 3 649 320 374];
CorrectX = input_data(x_input, S_steps, T_steps, L, T);

%% Iteration Scheme
xu(s(:), t(:)) = zeros(size(s,2), size(t,2));
xu = Boundary(CorrectX(s,1), CorrectU);

CorrectU(s(:), t(:)) = zeros(size(s,2), size(t,2));
xu2(s(:),t(:)) = zeros(size(s,2), size(t,2));
xu2 = Boundary(CorrectX(s,1), CorrectU);
    mkdir(folder_to_save);
    plotGraph2(xu2(1, :), xu(1, :), {0:T}, 't', 'x(s=0, t)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(2*S_steps/L+1, :), xu(2*S_steps/L+1, :), {0:T}, 't', 'x(s=2, t)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(5*S_steps/L+1, :), xu(5*S_steps/L+1, :), {0:T}, 't', 'x(s=5, t)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(end, :), xu(end, :), {0:T}, 't', 'x(s=7, t)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(:,1), xu(:,1), {0:L}, 's', 'x(s, t=0)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(:,2*T_steps/T+1), xu(:, 2*T_steps/T+1), {0:L}, 's', 'x(s, t=2)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(:,5*T_steps/T+1), xu(:, 5*T_steps/T+1), {0:L}, 's', 'x(s, t=5)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(:,end), xu(:,end), {0:L}, 's', 'x(s, t=8)', folder_to_save, 'Zero control', 'Non zero control');

%{



%% Init last column
%need_x(s) = xu(s, end);

need_x(s(1:(S_steps/L):end)) = [813983 150629 52499 20628 9938 3524 986 374]';
need_x(s) = interp1(s(1:(S_steps/L):end), need_x(s(1:(S_steps/L):end)), s(1:end));

%% Init of functionals
J2uCorrect = trapz(0:Ds:L, (xu(s, end) - need_x(s)').^2);
J1uCorrect = -p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectU(s, t)));

%% Optimization problem
pause
figure(1)
hold on;
xlabel('J2u')
ylabel('J1u')
lambda1 = 0;
while (lambda1 <= 1)
    k = 1;
    J2u = 1;
    uk = 5000*ones(size(s,2), size(t,2));
    uk_ = zeros(size(s,2), size(t,2));
    uprev = -1*ones(size(s,2), size(t,2));
    prev = 1;
    storedJ1u = [];
    storedJ2u = [];
    storedPrev = [];
    
    while (k < 1000000) && (prev > 10e-3) && (J2u > 10e-3)
        %% First Step
        % Pryamaya zadacha
        xu = Boundary(CorrectX(s,1), uk);

        % Obratnaya
        psiStartCondition = -2*(xu(s, end) - need_x(s)');
        psi = reverseBoundary(psiStartCondition);

        % Derivatives of Functionals
        J1_ = -p*phi*L*T;
        J2_ = psi;

        beta = sqrt(Ds^2+Dt^2)/norm(J2_);
        if beta < 0.1
            beta = 0.1;
        end

        uk_ = projection(uk - beta * (lambda1 * J1_ + lambda2 * J2_));

        %% Second step
        % Pryamaya zadacha
        xu = Boundary(CorrectX(s,1), uk_);

        % Obratnaya
        psiStartCondition = -2*(xu(s, end) - need_x(s)');
        psi = reverseBoundary(psiStartCondition);

        % Derivatives of Functionals
        J1_ = -p*phi*L*T;
        J2_ = psi;

        beta = k*(Ds^2+Dt^2)/norm(J2_);
        if beta < 0.1
            beta = 0.1;
        end

        uk = projection(uk_ - beta * (lambda1 * J1_ + lambda2 * J2_));


        xu = Boundary(CorrectX(s, 1), uk);

        k

        J2u = trapz((xu(s, end) - need_x(s)').^2)*Ds;
        J1u = -p*trapz(trapz(phi*uk(s, t))*Ds)*Dt;
        abs(J2u - J2uCorrect)
        abs(J1u - J1uCorrect)
        storedJ1u(end+1) = J1u;
        storedJ2u(end+1) = J2u;
        k = k+1;

        prev = sum(sum((uk-uprev).^2))
        storedPrev(end+1) = prev;
        
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
    folder_to_save = num2str(lambda1);
    mkdir(folder_to_save);
    %{
    plotGraph(CorrectX(1,:),xu(1,:), {0:T}, 't', 'x(s=0, t)', folder_to_save);
    plotGraph(CorrectX(1*S_steps/L+1,:),xu(1*S_steps/L+1,:), {0:T}, 't', 'x(s=1, t)', folder_to_save);
    plotGraph(CorrectX(2*S_steps/L+1,:),xu(2*S_steps/L+1,:), {0:T}, 't', 'x(s=2, t)', folder_to_save);
    plotGraph(CorrectX(5*S_steps/L+1,:),xu(5*S_steps/L+1,:), {0:T}, 't', 'x(s=5, t)', folder_to_save);
    plotGraph(CorrectX(end,:),xu(end,:), {0:T}, 't', 'x(s=7, t)', folder_to_save);
    
    plotGraph(CorrectX(:,1),xu(:,1), {0:L}, 's', 'x(s, t=0)', folder_to_save);
    plotGraph(CorrectX(:,1*T_steps/T+1),xu(:,1*T_steps/T+1), {0:L}, 's', 'x(s, t=1)', folder_to_save);
    plotGraph(CorrectX(:,2*T_steps/T+1),xu(:,2*T_steps/T+1), {0:L}, 's', 'x(s, t=2)', folder_to_save);
    plotGraph(CorrectX(:,5*T_steps/T+1),xu(:,5*T_steps/T+1), {0:L}, 's', 'x(s, t=5)', folder_to_save);
    plotGraph(CorrectX(:,t(end)),xu(:,t(end)), {0:L}, 's', 'x(s, t=8)', folder_to_save);

    plotGraph(CorrectU(1,:),uk(1,:), {0:T}, 't', 'u(s=0, t)', folder_to_save);
    plotGraph(CorrectU(1*S_steps/L+1,:),uk(1*S_steps/L+1,:), {0:T}, 't', 'u(s=1, t)', folder_to_save);
    plotGraph(CorrectU(2*S_steps/L+1,:),uk(2*S_steps/L+1,:), {0:T}, 't', 'u(s=2, t)', folder_to_save);
    plotGraph(CorrectU(5*S_steps/L+1,:),uk(5*S_steps/L+1,:), {0:T}, 't', 'u(s=5, t)', folder_to_save);
    plotGraph(CorrectU(end,:),uk(end,:), {0:T}, 't', 'u(s=7, t)', folder_to_save);

    plotGraph(CorrectU(:,1),uk(:,1), {0:L}, 's', 'u(s, t=0)', folder_to_save);
    plotGraph(CorrectU(:,1*T_steps/T+1),uk(:,1*T_steps/T+1), {0:L}, 's', 'u(s, t=1)', folder_to_save);
    plotGraph(CorrectU(:,2*T_steps/T+1),uk(:,2*T_steps/T+1), {0:L}, 's', 'u(s, t=2)', folder_to_save);
    plotGraph(CorrectU(:,5*T_steps/T+1),uk(:,5*T_steps/T+1), {0:L}, 's', 'u(s, t=5)', folder_to_save);
    plotGraph(CorrectU(:,t(end)),uk(:,t(end)), {0:L}, 's', 'u(s, t=8)', folder_to_save);
    save([pwd '/' folder_to_save '/storedJ1u.mat'], 'storedJ1u');
    save([pwd '/' folder_to_save '/storedJ2u.mat'], 'storedJ2u');
    save([pwd '/' folder_to_save '/storedPrev.mat'], 'storedPrev');
    %}
    lambda1 = lambda1+0.1;
    lambda2 = 1 - lambda1;
end
hold off;
saveas(gcf,'J1J2.jpg','jpg')

%}