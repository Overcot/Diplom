function [table] = input_data(data, S_steps, T_steps, L, T)
    global s t
    table(s(:), t(:)) = zeros(size(s,2), size(t,2));

    table(s(1:(S_steps/L):end), t(1:(T_steps/T):end)) = data;
    for j=t
        table(s, j) = interp1(s(1:(S_steps/L):end),table(s(1:(S_steps/L):end), j),s(1:end));
    end
    for i = s
        table(i, t) = interp1(t(1:(T_steps/T):end),table(i, t(1:(T_steps/T):end)),t(1:end));
    end
end