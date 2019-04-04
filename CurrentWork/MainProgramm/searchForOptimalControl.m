function [xOptim, uOptim, JOptim] = searchForOptimalControl(xData, uData, x0Data, L, T)
    %% Optimal Control Parameters
    global Ds Dt s t;

    rho = 0.3;
    p = 1;

    eps = 0;
    uMax = 1;
    uMin = 0.3;

    JData = trapz(0:Dt:T, trapz(0:Ds:L, functionalInternal(rho, p, xData, uData)))
    % for index = uMin:0.1:uMax
        index = uMin;
        uK = 1 * ones(size(s,2),size(t,2));
        k = 0;
        prev = 1;
        storedJu = [];
        storedPrev = [];
        uPrev = -1 * ones(size(s,2),size(t,2));

        while (k < 10) %&& (prev > 10e-3)
            %% First Step
            % Pryamaya zadacha
            xU = Boundary(x0Data, uK);
    
            % Obratnaya
            psi = reverseBoundary(xU, uK, L, rho);
    
            % Derivatives of Functionals
            J_ = JDerivative(xU, psi, rho, p);
    
            beta = calculateStep(k, 1);
    
            uK_ = projectionU(uK - beta * J_, uMin, uMax);
    
            %% Second step
            % Pryamaya zadacha
            xU = Boundary(x0Data, uK_);
    
            % Obratnaya
            psi = reverseBoundary(xU, uK, L, rho);
    
            % Derivatives of Functionals
            J_ = JDerivative(xU, psi, rho, p);
    
            beta = calculateStep(k, 1);
    
            uK = projectionU(uK_ - beta * J_, uMin, uMax);

            xU = Boundary(x0Data, uK);
    
            k
    
            % Ju = trapz(0:Dt:T, trapz(0:Ds:L, exp(-rho*t).*p.*uK(s, t).*xU(s, t)))
            Ju = sum(sum(exp(-rho*t).*p.*uK(s, t).*xU(s, t)))
            abs(Ju - JData)
            storedJu(end+1) = Ju;
            k = k+1;
    
            prev = sum(sum((uK-uPrev).^2))
            storedPrev(end+1) = prev;
            
            uPrev = uK;
        end
        uOptim = uPrev;
        xOptim = xU;
        JOptim = Ju;
        storedJu
        storedPrev
    % end
end

function [ J_ ] = JDerivative(x, psi, rho, p)
    global Ds Dt s t;
    global mu;
    % J' = -(1-mu) x * psi + delta'_u
    J_ = - (1 - mu(s))' .* x(s, t) .* psi(s, t) + functionalInternalDerivativeU(rho, p, x)
end

function [ psi ] = reverseBoundary(xU, uK, L, rho)
    global Ds Dt s t;
    global mu gamma;
    psi = zeros(size(s,2), size(t,2));
    right = M(1, uK, rho);
    for i=2:L
        right = right .* (1 + N(i, uK)) + M(i, uK, rho);
    end
    P = phiDerivative(gamma*xU);
    K = sum(gamma);
    A = K * P
    left = (1 + N(1, uK) - A);
    for i=2:L-1
        left = (1 + N(i, uK)).*left - A;
    end
    left = A - (1 + N(L, uK)).*left;

    psi0 = left.\right;
    psi(1, t) = psi0;
    int = sum(K*psi0);
    for class=s(end-1:-1:2)
        for time=t(class:end)
            psi(class, time) = MTime(class - 1, uK, rho, time - 1) + (1 + NTime(class, uK, time)).*psi(class - 1, time - 1) - int;
        end
    end
    psi
end

function value = M(i, u, rho)
    global s t;
    value = -exp(-rho*t).*u(i, t)
end

function value = MTime(i, u, rho, time)
    value = -exp(-rho*time).*u(i, time)
end

function value = N(i, u)
    global s t;
    global mu;
    value = mu(i)+(1 - mu(i)).*u(i, t);
end
function value = NTime(i, u, time)
    global mu;
    value = mu(i)+(1 - mu(i)).*u(i, time)
end
function [ nu_x ] = rightSideDerivativeX(u)
    global s t;
    global mu;
    nu_x = - mu(s) - (1 - mu(s)).*u(s, t)
end

function [ nu_u ] = rightSideDerivativeU(x)
    global s t;
    global mu;
    nu_u = - (1 - mu(s)) .* x(s, t)
end

function [ Jmatrix ] = functionalInternal(rho, p, xData, uData)
    global s t;
    Jmatrix = exp(-rho*t).*p.*uData(s, t).*xData(s, t);
end

function [ ans ] = functionalInternalDerivativeX(rho, p, uData)
    global s t;
    ans = exp(-rho*t).*p.*uData(s, t);
end

function [ ans ] = functionalInternalDerivativeU(rho, p, xData)
    global s t;
    ans = exp(-rho*t).*p.*xData(s, t);
end

function [ phi_ ] = phiDerivative(ksi)
    global s t;
    global a b allee ssbMax;
    phi_ = -1;
    if (ksi./ssbMax > allee)
        phi_ = exp(a-b.*ksi).*(1-b.*ksi);
    else
        phi_ = exp(ksi./ssbMax - allee) .* exp(a-b.*ksi) .* (1 - b.*ksi - ksi./ssbMax);
    end
end

function [ u ] = projectionU(u, uMin, uMax)
    global s t;
    for time=t
         for class=s
             if u(class, time) > uMax
                 u(class, time) = uMax;
             elseif u(class, time) < uMin
                 u(class, time) = uMin;
             end
         end
    end
 end

 function beta = calculateStep(k, JDerivative)
     global Ds Dt;
    %  beta = k * (Ds^2 + Dt^2)/norm(JDerivative);
    %  if beta < 0.1
    %      beta = 0.1;
    %  end
    beta = 0.1;
 end
 