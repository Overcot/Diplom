clear; clc;
global recruitment ssb

recruitment = xlsread('SRdata_ICES','D:D')
ssb = xlsread('SRdata_ICES','C:C')

init = [1, 1];

options = optimset('Display','iter','PlotFcns',@optimplotfval,'Tolx',1e-20,'TolFun',1e-20);
[x, val, exitflag,output] = fminsearch(@recruitment_func, init,options)

x(1)
x(2)