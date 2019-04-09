function [xOptim, uOptim, JOptim, storedJu] = searchForOptimalControl(xData, uData, x0Data, L, T)
    %% Optimal Control Parameters
    global Ds Dt s t;

    rho = 0.3;
    p = 1;

    eps = 0;
    uMax = 1;
    uMin = 0.35;

    % JData = sum(sum(functionalInternal(rho, p, xData, uData)))
    % for index = uMin:0.1:uMax
        index = uMax;
        uK = (uMin+uMax)/2 * ones(size(s,2),size(t,2));
        k = 0;
        prev = 1;
        storedJu = [];
        storedPrev = [];
        uPrev = -1 * ones(size(s,2),size(t,2));

        while (k < 10000) %&& (prev > 10e-3)
            %% First Step
            % Pryamaya zadacha
            xU = Boundary(x0Data, uK);
    
            % Obratnaya
            psi = reverseBoundary(xU, uK, L, T, rho);
            
            % Derivatives of Functionals
            J_ = JDerivative(xU, psi, rho, p);
    
            beta = calculateStep(k, J_);
    
            uK_ = projectionU(uK + beta * J_, uMin, uMax);
    
            %% Second step
            % Pryamaya zadacha
            xU_ = Boundary(x0Data, uK_);
    
            % Obratnaya
            psi_ = reverseBoundary(xU_, uK_, L, T, rho);
    
            % Derivatives of Functionals
            J__ = JDerivative(xU_, psi_, rho, p);
    
            beta_ = calculateStep(k, J__);
    
            uK = projectionU(uK + beta_ * J__, uMin, uMax);

            xU = Boundary(x0Data, uK);
    
            k
    
            % Ju = trapz(0:Dt:T, trapz(0:Ds:L, exp(-rho*t).*p.*uK(s, t).*xU(s, t)))
            Ju = sum(sum(functionalInternal(rho, p, xU, uK)))
            storedJu(end+1) = Ju;
            k = k+1;
    
            prev = sum(sum((uK-uPrev).^2));
            storedPrev(end+1) = prev;
            
            uPrev = uK;
        end
        uOptim = uPrev;
        xOptim = xU;
        JOptim = Ju;
        storedJu;
        storedPrev;
    % end
end

function [ J_ ] = JDerivative(x, psi, rho, p)
    global Ds Dt s t;
    global mu;
    % J' = -(1-mu) x * psi + delta'_u
    J_ = - (1 - mu(s))' .* x(s, t) .* psi(s, t) + functionalInternalDerivativeU(rho, p, x);
end

function [ psi ] = reverseBoundary(xU, uK, L, T, rho)
    global Ds Dt s t;
    global mu gamma;
    psi = zeros(size(s,2), size(t,2));
    left = zeros((size(s,2) - 1)*(size(t,2) - 1));
    right = zeros((size(s,2) - 1)*(size(t,2) - 1), 1);

    K = sum(gamma);
    P = phiDerivative(gamma*xU);
    KPHI = K * P(1:end-1);
    for time=t(1:end-1)
        left(:,time) = KPHI(time);
    end

    for class=s(1:end-1)
        for time=t(1:end-1)
            i = (class-1)*(size(t,2)-1) + time;
            left(i, i) = left(i, i) - (1+NTime(class, uK, time));

            right(i) = MTime(class, uK, rho, time);
        end
    end

    matrixOnes = diag(ones(size(t,2)-2, 1), 1);
    for index=1:L-1 
        left((index-1)*T+1:index*T,index*T+1:(index+1)*T) = matrixOnes;
    end
    
    X = mldivide(left, right);

    for class=s(1:end-1)
        for time=t(1:end-1)
            (class-1)*(size(t,2)-1) + time;
            psi(class, time)=X((class-1)*(size(t,2) - 1) + time);
        end
    end
end

function value = M(i, u, rho)
    global s t;
    value = -exp(-rho*(t-1)).*u(i, t);
end

function value = MTime(i, u, rho, time)
    value = -exp(-rho*(time-1)).*u(i, time);
end

function value = N(i, u)
    global s t;
    global mu;
    value = mu(i)+(1 - mu(i)).*u(i, t);
end

function value = NTime(i, u, time)
    global mu;
    value = mu(i)+(1 - mu(i)).*u(i, time);
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

function [ Jmatrix ] = functionalInternal(rho, p, x, u)
    global s t;
    Jmatrix = exp(-rho*(t-1)).*p.*u(s, t).*x(s, t);
end

function [ ans ] = functionalInternalDerivativeX(rho, p, u)
    global s t;
    ans = exp(-rho*(t-1)).*p.*u(s, t);
end

function [ ans ] = functionalInternalDerivativeU(rho, p, x)
    global s t;
    ans = exp(-rho*(t-1)).*p.*x(s, t);
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
    beta = 0.01*1/norm(JDerivative);
 end
 