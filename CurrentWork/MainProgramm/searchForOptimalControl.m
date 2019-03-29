function [xOptim, uOptim] = searchForOptimalControl(xData, uData, x0Data, L, T)
    %% Optimal Control Parameters
    global Ds Dt s t;

    rho = 0.3;
    p = 1;

    eps = 0;
    uMax = 1;
    uMin = 0;

    JData = trapz(0:Dt:T, trapz(0:Ds:L, exp(-rho*t).*p.*uData(s, t).*xData(s, t)))
    for index = uMin:0.1:uMax

        uK = ones(size(s,2),size(t,2));

        while (k < 10) && (prev > 10e-3) && (J2u > 10e-3)
            %% First Step
            % Pryamaya zadacha
            xu = Boundary(x0Data, uK);
    
            % Obratnaya
            psiStartCondition = -2*(xu(s, end) - need_x(s)');
            psi = reverseBoundary(psiStartCondition);
    
            % Derivatives of Functionals
            
            % J_ = psi;
    
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
    end
    
    
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
    %{
    lambda1 = lambda1+0.1;
    lambda2 = 1 - lambda1;
    %}

end