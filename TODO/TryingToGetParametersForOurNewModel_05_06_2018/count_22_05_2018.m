clear; clc;
global recruitment ssb

%% Import Data
year = 2011
if (year == 2010)
    recruitment = xlsread('SRdata_ICES','D2:D49')
    ssb = xlsread('SRdata_ICES','C2:C49')
elseif (year == 2011)
    recruitment = xlsread('SRdata_ICES','D2:D50')
    ssb = xlsread('SRdata_ICES','C2:C50')
end

init = [1000, 0];

options = optimset('Display','iter','PlotFcns',@optimplotfval,'Tolx',1e-20,'TolFun',1e-20);
[x, val, exitflag, output] = fminsearch(@recruitment_func, init,options);
x(1)
x(2)

%% Curve Fit
xdata = ssb;
ydata = log(recruitment);
ricker(init, xdata);
%%[x, resnorm] = lsqcurvefit(@ricker, init, xdata, ydata)

%% Count AIC
resnorm = 11.44
n = size(recruitment, 1)
m = 2;
L = -(n/2.0)*(log(2*pi)+log(resnorm/n)+1);
AIC = -2*L+2*m