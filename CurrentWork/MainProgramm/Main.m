clear all; 
clc;
global ssb;
global x Ds Dt s t alpha 
global gamma mu;
global a b allee;
global umin umax;
precision = 1;
L = 10;
T = 50;
gammaValue = [0.013622642 0.161377358 0.419283019 0.698283019 0.877773585 0.999245283 1 1 1 1 1];
muValue = [0.319 0.409 0.788 0.818 0.818 0.818 0.818 0.818 0.818 0.818 0.818];
x0 = [487478000 184241000 20339000 10330000 8251000 871500 871500 871500 871500 871500 871500];
%anna model
a = 8.244595391;
b = 8.58578E-07;
allee = 0.2;

folder_to_save = ['graphs,precision=(',num2str(precision), ')'];


%% Init Part
[x, u, Ds, Dt, S_steps, T_steps, s, t, gamma, mu] = startInit(precision, L, T, gammaValue, muValue, x0);

x = Boundary(x0,u);

%% Params of Setki

lambda1 = 0;
lambda2 = 1 - lambda1;
phi = 1;
p = 1;

eps = 0;
umax = 10000;
umin = 0;

mkdir(folder_to_save);
plotGraph(x(1,:), {0:T}, 't', 'x(s=0,t)', 'Zero control', folder_to_save);
%{
    mkdir(folder_to_save);
    plotGraph2(xu2(1, :), xu(1, :), {0:T}, 't', 'x(s=0, t)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(2*S_steps/L+1, :), xu(2*S_steps/L+1, :), {0:T}, 't', 'x(s=2, t)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(5*S_steps/L+1, :), xu(5*S_steps/L+1, :), {0:T}, 't', 'x(s=5, t)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(end, :), xu(end, :), {0:T}, 't', 'x(s=7, t)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(:,1), xu(:,1), {0:L}, 's', 'x(s, t=0)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(:,2*T_steps/T+1), xu(:, 2*T_steps/T+1), {0:L}, 's', 'x(s, t=2)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(:,5*T_steps/T+1), xu(:, 5*T_steps/T+1), {0:L}, 's', 'x(s, t=5)', folder_to_save, 'Zero control', 'Non zero control');
    plotGraph2(xu2(:,end), xu(:,end), {0:L}, 's', 'x(s, t=8)', folder_to_save, 'Zero control', 'Non zero control');
%}





%% Init of functionals
%J2uCorrect = trapz(0:Ds:L, (xu(s, end) - need_x(s)').^2);
%J1uCorrect = -p*trapz(0:Dt:T, trapz(0:Ds:L, phi*CorrectU(s, t)));

%% Optimization problem
%{
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
