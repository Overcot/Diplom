function x0 = recruitmentFunction(model, time, x, params)
%   function to count x(0,t) for any given t
%       input:
%           model - suitable options: 'Tahvonen' and 'New'(to be implemented)
%           time - we search for x(1, time)
%           x - all matrix of solution
%           params - [Ds, Dt, (params for each model)]
%               for Tahvonen: (gamma, p(time))
%               for New: ...
%       output:
%           x(0,time) = -1 (in matlab its x(1,time))(-1 by default for
%           error check(to be implemented)
%           
%
    x0 = -1;
    Ds = params(1); Dt = params(2);
    if model == 'Tahvonen'
        gamma = params(3); p = params(4)
        x0 = (p + trapz(gamma*x(2:end, time)*Ds))/(1-gamma*Ds);
    elseif model == 'New'
        return
    end
end