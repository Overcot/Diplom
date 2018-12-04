function x0 = recruitmentFunction(model, time, x, params)
global ssb
global a b allee gamma p
%   function to count x(0,t+1) for any given t
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
    if strcmp(model,'Tahvonen')
        %x0 = (p + trapz(gamma*x(2:end, time)*Ds))/(1-gamma*Ds);
    elseif strcmp(model,'Anna')
        ssbCurr = gamma*x(1:end, time)*Ds/1000;
        
        ssb = [ssb, ssbCurr];
        if ssb(end-1)./274032 > allee
            value = log(ssb(end-1)) + a - b*ssb(end-1)
        else
            value = log(ssb(end-1)) + a - b*ssb(end-1) + (ssb(end-1)./274032-allee)
        end
        x0 = exp(value);
    end
end