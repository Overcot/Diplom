function [ psi ] = reverseBoundary( startcondition )
    global Ds Dt s t mu f gamma;
    psi = zeros(size(s,2), size(t,2));
    psi(:, end) = startcondition;
    for time=t(end:-1:2)
        psi(s(1:end - 1) + 1, time - 1)
        (Ds/Dt)*(psi(s(1:end - 1) + 1, time) - psi(s(1:end - 1) + 1, time - 1))
        (mu(s(1:end-1)+1)-1)'.*psi(s(1:end - 1) + 1, time) + gamma(s(1:end - 1) + 1).* psi(1, time)
        
        psi(s(1:end - 1), time - 1) = psi(s(1:end - 1) + 1, time - 1) + (Ds/Dt)*(psi(s(1:end - 1) + 1, time) - psi(s(1:end - 1) + 1, time - 1)) + Ds*( (mu(s(1:end-1)+1)-1)'.*psi(s(1:end - 1) + 1, time) + gamma(s(1:end - 1) + 1) .* psi(1, time));
    end
end

