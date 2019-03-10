function x0 = recruitmentFunction(model, time, x, params)
global ssb
global a b allee gamma
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
        x0 = (gamma(2:end)*x(2:end, time)*Ds)/(1-gamma(1)*Ds);
    elseif strcmp(model,'Anna')
        %{
        ssbCurr = gamma*x(1:end, time)*Ds;        
        ssb= [ssb, ssbCurr];
        if ssbCurr./274032 > allee
            value = log(ssbCurr) + a - b*ssbCurr
        else
            value = log(ssbCurr) + a - b*ssbCurr + (ssbCurr./274032-allee)
        end
        x0 = exp(value);
        %}
        I = gamma(2:end)*x(2:end, time)*Ds;
        x0 = exp(a-b*I)*I/(1-gamma(1)*Ds*exp(a-b*I))
    end
end