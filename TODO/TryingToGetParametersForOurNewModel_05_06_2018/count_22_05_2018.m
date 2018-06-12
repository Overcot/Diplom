clear; clc;
global recruitment ssb

recruitment = xlsread('SRdata_ICES','D3:D50')
ssb = xlsread('SRdata_ICES','C2:C49')

init = [3416, 5*10^(-8)];

options = optimset('Display','iter','PlotFcns',@optimplotfval,'Tolx',1e-20,'TolFun',1e-20);
[x, val, exitflag,output] = fminsearch(@recruitment_func, init,options)

x(1)
x(2)