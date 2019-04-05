function x0 = recruitmentFunction(model, time, x)
global ssbMax
global a b allee gamma
global Ds Dt
%   function to count x(0, t + 1) for any given t
%       input:
%           model - suitable options: 'Tahvonen' and 'Anna'
%           time - we search for x(1, time)
%           x - all matrix of solution
%       output:
%           x(0,time) = -1 (in matlab its x(1,time))(-1 by default for
%           error check(to be implemented)
%           
%
    x0 = -1;
    if strcmp(model, 'Tahvonen')
        x0 = (gamma(2:end)*x(2:end, time)*Ds)/(1-gamma(1)*Ds);
    elseif strcmp(model, 'Anna')
        ssbCurr = gamma(1:end)*x(1:end, time)*Ds;
        if (ssbCurr./ssbMax > allee)
            x0 = ssbCurr*exp(a-b*ssbCurr);
        else
            x0 = ssbCurr*exp(a-b*ssbCurr)*exp(ssbCurr./ssbMax - allee);
        end
    end
end