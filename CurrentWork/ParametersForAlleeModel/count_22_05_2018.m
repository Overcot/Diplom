clear; clc;
global recruitment ssb maxSSB lowerSSB higherSSB

year = 2011 %possible values: 2011, 2010
model = 'Anna' %possible values: 'Anna', 'Ricker'

%% Import Data
if (year == 2010)
    recruitment = xlsread('ourData2010&2011','D3:D49');
    ssb = xlsread('ourData2010&2011','B4:B50');
    maxSSB = max(ssb);
elseif (year == 2011)
    recruitment = xlsread('ourData2010&2011','H3:H50');
    ssb = xlsread('ourData2010&2011','F4:F51');
    maxSSB = max(ssb);
end

options = optimset('Display','iter','PlotFcns',@optimplotfval,'Tolx',1e-20,'TolFun',1e-20);
%% Try finding parameters using fminsearch
if (model == 'Anna')
    init = [0, 0, 0];
    lowerSSB = 0;
    higherSSB = 0;
    [x, val, exitflag, output] = fminsearch(@AnnaModel, init, options);
elseif (model == 'Ricker')
    init = [0, 0];
    lowerSSB = 0;
    higherSSB = 0;
    [x, val, exitflag, output] = fminsearch(@recruitment_func, init, options);
end
x
val
lowerSSB
higherSSB
%% Try finding parameters using fmincon
% Init data
xdata = ssb;
ydata = log(recruitment);
optionslsqnonlin = optimoptions(@lsqnonlin, 'Algorithm','trust-region-reflective', 'Display', 'iter', 'Tolx',1e-20,'TolFun',1e-20);
A = [];
b = [];
Aeq = [];
beq = [];
if (model == 'Anna')
    init = [8,0,0.14];
    lowerSSB = 0;
    higherSSB = 0;
    [x, resnorm] = fmincon(@AnnaModel, init, A,b,Aeq,beq,[-10000000,-10000000,0], [100000000,10000000,0.25]);
elseif (model == 'Ricker')
    init = [3000, 0];
    lowerSSB = 0;
    higherSSB = 0;
    [x, resnorm] = lsqnonlin(@ricker, init,[],[], options);
end

x
resnorm
lowerSSB
higherSSB
fileName = 'possibleVariants.xlsx'
N='Adnan'; a=22; roll=22; gpa=3.55;
fileExist = exist(fileName,'file'); 
if fileExist==0
    header = {'Name', 'age ','roll' , 'gpa'};
    xlswrite(filename,header);


else

    [~,~,input] = xlsread(filename); % Read in your xls file to a cell array (input)
    new_data = {N, a,roll , gpa}; % This is a cell array of the new line you want to add
    output = cat(1,input,new_data); % Concatinate your new data to the bottom of input
    xlswrite(filename,output); % Write to the new excel file. 

end
%% Count AIC
%resnorm = 11.44
%n = size(recruitment, 1)
%m = 2;
%L = -(n/2.0)*(log(2*pi)+log(resnorm/n)+1);
%AIC = -2*L+2*m