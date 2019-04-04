clear all; clc;
global ssbMax;
global x Ds Dt s t;
global gamma mu;
global a b allee;
global umin umax;

%% Установка констант
excelToImportFrom = 'AllDataICES2016'
precision = 1
L = 5
T = 20
baseStartingYear = 1963;
ourStartingYear = 1964

%% Import data from excel
[xData, x0Data, ssbData, ssbMax, fishMortalityData, gammaData, muData, a, b, allee] = importFromExcel(excelToImportFrom, ourStartingYear, baseStartingYear, T);

%% Init Part based on constants and data
[x, x0, u, Ds, Dt, S_steps, T_steps, s, t, gamma, mu] = startInit(precision, L, T, gammaData, muData, x0Data, fishMortalityData);

x = Boundary(x0, u);

allee = 0
xAlleeZero = Boundary(x0, u);



% folder_to_save = 'WithAndWithoutAllee';
% if ~exist(folder_to_save, 'dir')
%     mkdir(folder_to_save)
% end
% for i=1:L+1
%     msg = strcat('x(s=',num2str(i),',t)');
%     plotGraph2(xAlleeZero(i,:), x(i,:), T+1, 't', msg, 'WithoutAllee', 'Allee', folder_to_save, 'Pop numbers (In Thousands) with Allee Vs Without Allee');
% end

% for i=1:T+1
%     msg = strcat('x(s,t=',num2str(i),')');
%     plotGraph2(xAlleeZero(:,i), x(:,i), L+1, 's', msg, 'WithoutAllee', 'Allee', folder_to_save, 'Pop numbers (In Thousands) with Allee Vs Without Allee');
% end

% folder_to_save = 'WithAllee_VS_Data';
% if ~exist(folder_to_save, 'dir')
%     mkdir(folder_to_save)
% end

% for i=1:L+1
%     msg = strcat('x(s=',num2str(i),',t)');
%     plotGraph2(xData(i,:), x(i,:), T+1, 't', msg, 'data', 'numerical', folder_to_save, 'Pop numbers (In Thousands) with Allee Vs Data');
% end

% for i=1:T+1
%     msg = strcat('x(s,t=',num2str(i),')');
%     plotGraph2(xData(:,i), x(:,i), L+1, 's', msg, 'data', 'numerical', folder_to_save, 'Pop numbers (In Thousands) with Allee Vs Data');
% end

[xOptim, uOptim, JOptim] = searchForOptimalControl(xData, fishMortalityData, x0Data, L, T);
rho = 0.3;
p = 1;
JData = sum(sum(exp(-rho*t).*p.*fishMortalityData(s, t).*xData(s, t)))
folder_to_save = 'OptimalSolution';
if ~exist(folder_to_save, 'dir')
    mkdir(folder_to_save)
end

for i=1:L+1
    msg = strcat('x(s=',num2str(i),',t)');
    plotGraph2(xData(i,:), xOptim(i,:), T+1, 't', msg, 'data', 'optimal', folder_to_save, 'Pop numbers In Optimal Solution (In Thousands) with Allee Vs Data');
end

for i=1:T+1
    msg = strcat('x(s,t=',num2str(i),')');
    plotGraph2(xData(:,i), xOptim(:,i), L+1, 's', msg, 'data', 'optimal', folder_to_save, 'Pop numbers In Optimal Solution (In Thousands) with Allee Vs Data');
end

for i=1:L+1
    msg = strcat('u(s=',num2str(i),',t)');
    plotGraph2(fishMortalityData(i,:), uOptim(i,:), T+1, 't', msg, 'data', 'optimal', folder_to_save, 'Fish Mortality with Allee Vs Data Fish Mortality');
end

for i=1:T+1
    msg = strcat('u(s,t=',num2str(i),')');
    plotGraph2(fishMortalityData(:,i), uOptim(:,i), L+1, 's', msg, 'data', 'optimal', folder_to_save, 'Fish Mortality with Allee Vs Data Fish Mortality');
end