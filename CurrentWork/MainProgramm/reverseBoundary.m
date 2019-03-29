function [ psi ] = reverseBoundary( startcondition )
    global Ds Dt s t;
    global mu gamma;
    psi = zeros(size(s,2), size(t,2));
    psi(:, end) = 0;
    psi(end, :) = 0;
    for time=t(end:-1:2)
        for class = s(1:end-2)
            % psi's(class, time) + psi't(class, time) = - exp(-rho * time) * p * u(class, time) - psi(class, time)

        end



        psi(s(1:end - 1) + 1, time - 1)
        (Ds/Dt)*(psi(s(1:end - 1) + 1, time) - psi(s(1:end - 1) + 1, time - 1))
        (mu(s(1:end-1)+1)-1)'.*psi(s(1:end - 1) + 1, time) + gamma(s(1:end - 1) + 1).* psi(1, time)
        
        
        
        
        
        psi(s(1:end - 1), time - 1) = psi(s(1:end - 1) + 1, time - 1) + (Ds/Dt)*(psi(s(1:end - 1) + 1, time) - psi(s(1:end - 1) + 1, time - 1)) + Ds*( (mu(s(1:end-1)+1)-1)'.*psi(s(1:end - 1) + 1, time) + gamma(s(1:end - 1) + 1) .* psi(1, time));
    end
end

