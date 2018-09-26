clear; clc;
global recruitment ssb

%% Import Data
year = 2011
if (year == 2010)
    recruitment = xlsread('ourData2010&2011','D3:D49');
    ssb = xlsread('ourData2010&2011','B4:B50');
elseif (year == 2011)
    recruitment = xlsread('ourData2010&2011','H3:H50');
    ssb = xlsread('ourData2010&2011','F4:F51');
end

init = [0, 0,0];

options = optimset('Display','iter','PlotFcns',@optimplotfval,'Tolx',1e-20,'TolFun',1e-20);
%[x, val, exitflag, output] = fminsearch(@recruitment_func, init,options);
[x, val, exitflag, output] = fminsearch(@AnnaModel, init,options);
x(1)
x(2)
x(3)
val

%% Curve Fit
xdata = ssb;
ydata = log(recruitment);
optionslsqnonlin = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
[x, resnorm] = lsqnonlin(@ricker, init,[],[], options)
[ricker(x), zeros(size(xdata,1),1)]

%% Count AIC
resnorm = 11.44
n = size(recruitment, 1)
m = 2;
L = -(n/2.0)*(log(2*pi)+log(resnorm/n)+1);
AIC = -2*L+2*m