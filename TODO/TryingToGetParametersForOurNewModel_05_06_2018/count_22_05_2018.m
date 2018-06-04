clear; clc;
global recruitment ssb

recruitment = xlsread('SRdata_ICES','D:D')
ssb = xlsread('SRdata_ICES','C:C')

init = [1, 1];

options = optimset('Display','iter','PlotFcns',@optimplotfval,'Tolx',1e-20,'TolFun',1e-20);
[x, val, exitflag,output] = fminsearch(@recruitment_func, init,options)

min = -1;
ansA = 0;
ansB = 0;
%{
for a=900000:0.01:10000000
    a
    for b=0:0.01:100
        tmp = sum(recruitment(t)-ssb(t).*exp(a-b.*ssb(t)).^2);
        if isfinite(tmp)
            if min > 0
                if tmp < min
                    min = tmp;
                    ansA = a;
                    ansB = b;
                end
            else
                min = tmp;
                ansA = a;
                ansB = b;
            end
        end
    end
end
ansA
ansB

%ansA = 995352
%ansB = 24
%ansA = -821748
%ansB = -2
%ansA = 9.0536967
%ansB = 20.93
%}