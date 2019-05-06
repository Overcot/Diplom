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

% allee = 0
% xAlleeZero = Boundary(x0, u);
% plotGraph(1-mu(:), size(mu, 2), '�������', '�����������','�����������', '');

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
x = Boundary(x0Data, fishMortalityData);

if exist('results.mat')
    load('results.mat','xOptim','uOptim','storedJ1u','storedJ2u','storedL','storedLambda','storedK','storedPrev');
else
    xOptim = [];
    uOptim = [];
    storedJ1u = [];
    storedJ2u = [];
    storedL = [];
    storedLambda = [];
    storedK = 1;
    storedPrev = [];
end
addPoints = 15000;
[xOptim, uOptim, J1Optim, J2Optim, storedJ1u, storedJ2u, storedL, storedLambda, storedK, storedPrev] = searchForOptimalControl(addPoints, xData, fishMortalityData, x0Data, L, T, xOptim, uOptim, storedJ1u, storedJ2u, storedL, storedLambda, storedK, storedPrev);
save('results.mat', 'xOptim','uOptim','storedJ1u','storedJ2u','storedL','storedLambda','storedK', 'storedPrev')
rho = 0.3;
p = 1;
x = Boundary(x0Data, fishMortalityData);
storedL(end)
J1Data = sum(sum(exp(-rho*(t-1)).*p.*fishMortalityData(s, t).*x(s, t)));
J2Data = sum(sum(fishMortalityData.^2));

folder_to_save = 'OptimalSolution';
if ~exist(folder_to_save, 'dir')
    mkdir(folder_to_save)
end
%{
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
%}
for i=1:T+1
    name = strcat('lambda(',num2str(i),')');
    plotGraph(storedLambda(i,:), size(storedLambda, 2), 'index', name, 'lambda', folder_to_save);
end

plotGraph(gamma*xOptim, size(t,2), 't', 'SSB(t)','SSB(t)', folder_to_save);
plotGraph(storedJ1u(:), size(storedJ1u, 2), 'index', 'J1u(index)','Functional', folder_to_save);
plotGraph(storedJ2u(:), size(storedJ2u, 2), 'index', 'J2u(index)','Functional', folder_to_save);
plotGraph(storedL(:), size(storedL, 2), 'index', 'L(u, lambda)', 'Свертка', folder_to_save);
