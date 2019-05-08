function [xOptim, uOptim, J1Optim, J2Optim, storedJ1u, storedJ2u, storedL, storedLambda, k, storedPrev] = searchForOptimalControl(add, xData, uData, x0Data, L, T, xOptim, uOptim, storedJ1u, storedJ2u, storedL, storedLambda, storedK, storedPrev)
    %% Optimal Control Parameters
    global Ds Dt s t;
    global ssbMax allee;
    rho = 0.3;
    p = 1;

    eps = 0;
    uMax = 0.8;
    uMin = 0.3;

    lambdaMin = 0;
    % for index = uMin:0.1:uMax
        index = uMax;
        if isempty(uOptim)
            uK = (uMin + uMax)/2 * ones(size(s,2),size(t,2));
            uK(:,1) = 0;
        else
            uK = uOptim;
        end
        
        if ~isempty(xOptim)
            xU = xOptim;
        end
        
        if isempty(storedLambda)
            lambdaK = zeros(1, size(t,2));
        else
            lambdaK = storedLambda(:,end)';
        end
        
        k = storedK;
        prev = 1;
        storedJ1u = storedJ1u;
        storedJ2u = storedJ2u;
        storedPrev = storedPrev;
        storedL = storedL;
        storedLambda = storedLambda;
        uPrev = -1 * ones(size(s,2),size(t,2));
        while (k < add + storedK) %&& (prev > 10e-3)
            %% First Step
            % Pryamaya zadacha
            xU = Boundary(x0Data, uK);
    
            % Obratnaya
            psi = reverseBoundary(xU, uK, L, T, rho, lambdaK);
            
            % Derivatives of Functionals
            J1_ = JDerivative(xU, psi, rho, p);

            lambdaDer = lambdaDerivative(xU);

            % Calculate steps
            beta = calculateStep(k, J1_);
            alpha = 10^-6;

            % Prognoznii step
            uK_ = projectionU(uK - beta * J1_, uMin, uMax);
            lambdaK_ = projectionLambda(lambdaK + alpha*lambdaDer, lambdaMin);

            %% Second step
            % Pryamaya zadacha
            xU_ = Boundary(x0Data, uK_);
    
            % Obratnaya
            psi_ = reverseBoundary(xU_, uK_, L, T, rho, lambdaK_);
    
            % Derivatives of Functionals
            J1__ = JDerivative(xU_, psi_, rho, p);

            lambdaDer_ = lambdaDerivative(xU_);

            beta_ = calculateStep(k, J1__);

            % Main step
            uK = projectionU(uK - beta_ * J1__, uMin, uMax);
            lambdaK = projectionLambda(lambdaK + alpha*lambdaDer_, lambdaMin);


            xU = Boundary(x0Data, uK);
            k
            
            gU = sum(sum(functionalInternal(rho, p, xU, uK)));
            J1u = sum(sum(-functionalInternal(rho, p, xU, uK)))
            J2u = sum(sum(uK.^2));

            storedJ1u(end+1) = J1u;
            storedJ2u(end+1) = J2u;
            
            k = k+1;
    
            prev = sum(sum((uK - uPrev).^2));
            storedPrev(end+1) = prev;
            lastL = gU + sum(lambdaK.*lambdaDerivative(xU))
            storedL(end+1) = lastL;
            
            storedLambda(:,end+1) = lambdaK';

            uPrev = uK;
        end
        uOptim = uPrev;
        xOptim = xU;
        J1Optim = storedJ1u(end);
        J2Optim = storedJ2u(end);

        storedJ1u;
        storedJ2u;
        storedPrev;
        storedL;
    % end
end

function [ J_ ] = JDerivative(x, psi, rho, p)
    global Ds Dt s t;
    global mu;
    % J' = -(1-mu) x * psi + delta'_u
    J_ = rightSideDerivativeU(x) .* psi(s, t) + functionalInternalDerivativeU(rho, p, x);
end

function [ lambdaDerivative ] = lambdaDerivative(x)
    global s t;
    global gamma ssbMax allee;
    lambdaDerivative = allee*ssbMax - gamma(s)*x(s, t);
end

function [ psi ] = reverseBoundary(xU, uK, L, T, rho, lambda)
    global Ds Dt s t;
    global mu gamma;
    psi = zeros(size(s,2), size(t,2));
    left = zeros((size(s,2) - 1)*(size(t,2) - 1));
    right = zeros((size(s,2) - 1)*(size(t,2) - 1), 1);

    K = sum(gamma);
    P = phiDerivative(gamma*xU);
    KPHI = K * P(1:end-1);
    for time=t(1:end-1)
        left(:,time) = - KPHI(time);
    end

    for class=s(1:end-1)
        for time=t(1:end-1)
            i = (class-1)*(size(t,2)-1) + time;
            left(i, i) = left(i, i) + 1;

            right(i) = MTime(class+1, uK, rho, time+1) - lambda(time+1).*gamma(class+1);
        end
    end
    
    for index=1:L-1 
        matrixOnes = diag(ones(size(t,2)-2, 1), 1);
        for time=t(1:end-2)
            matrixOnes(time, time+1) = NTime(index, uK, time);
        end
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
    value = mu(i)+(1 - mu(i)).*u(i, t) - 1;
end

function value = NTime(i, u, time)
    global mu;
    value = mu(i)+(1 - mu(i)).*u(i, time) - 1;
end

function [ nu_x ] = rightSideDerivativeX(u)
    global s t;
    global mu;
    nu_x = - mu(s) - (1 - mu(s)).*u(s, t);
end

function [ nu_u ] = rightSideDerivativeU(x)
    global s t;
    global mu;
    nu_u = - (1 - mu(s)) * x(s, t);
end

function [ Jmatrix ] = functionalInternal(rho, p, x, u)
    global s t;
    Jmatrix = -exp(-rho*(t-1)).*p.*u(s, t).*x(s, t);
end

function [ ans ] = functionalInternalDerivativeX(rho, p, u)
    global s t;
    ans = -exp(-rho*(t-1)).*p.*u(s, t);
end

function [ ans ] = functionalInternalDerivativeU(rho, p, x)
    global s t;
    ans = -exp(-rho*(t-1)).*p.*x(s, t);
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

function lambda = projectionLambda(lambda, lambdaMin)
    global t
    for time=t
        if lambda(time) < lambdaMin
            lambda(time) = lambdaMin;
        end
    end
end

 function beta = calculateStep(k, JDerivative)
    % beta = 10*1/norm(JDerivative);
    beta = 10^-9;
 end
 