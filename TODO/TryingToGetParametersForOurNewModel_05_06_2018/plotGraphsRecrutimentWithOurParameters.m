clear; clc;
global recruitment ssb maxSSB lowerSSB higherSSB

year = 2011 %possible values: 2011, 2010
model = 'Anna' %possible values here: 'Anna'

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
tYears = 1963:year-1;

%% Params for Ricker with Allee

params1 = [8.244595391; 8.58578E-07; 0.2];
params2 = [8.635711156;	2.64626E-06; 0.54595332];

J1 = AnnaModel(params1)
J2 = AnnaModel(params2)

%% Plot graphs Rec(t)

y = log(recruitment);
y1 = AnnaModelRecruitment(params1);
y2 = AnnaModelRecruitment(params2);

figure;
hold on;
plot(tYears, y, '*-');
plot(tYears, y1, tYears, y2);
legend('data/recruitment', 'a = 8.24459539, b = 8.58578E-07, allee = 0.2','a = 8.635711156, b = 2.64626E-06, allee = 0.54595332');
xlabel('Time(Years)')
ylabel('ln (recruitment)')
title('Comparison ICES 2011 & Ricker With Allee Effect Function with different parametrs');

%% graph of allee effect - (1/x)*(Dx/Dt) of (x)
t=1:length(recruitment)-1; 

% graph of data 
[yDataSorted,indexOfYDataSorted] = sort(y(t),'ascend');
dydt = (y(t+1)-y(t))./y(t);
dydtSorted = dydt(indexOfYDataSorted);

% graph of recruitment with allee
[y1DataSorted, indexOfY1DataSorted] = sort(y1(t), 'ascend');
dy1dt = (y1(t+1)-y1(t))./y1(t);
dy1dtSorted = dy1dt(indexOfY1DataSorted);

%{
figure;
hold on;
plot(yDataSorted, dydtSorted);
plot(y1DataSorted, dy1dtSorted);

xlabel('R');
ylabel('1/R * dR/dt');
legend('data','a = 8.24459539, b = 8.58578E-07, allee = 0.2');
%}

%% Graph of R/SSb(Ssb) - similar in PhD_seminar_Allee_effect.pdf(Page 19)

[ssbSorted, indexOfSSbSorted] = sort(ssb(t), 'ascend'); % sorting ssb and saving indexes

ssbSorted1 = ssbSorted*100/(max(ssbSorted)); % comment to show in real values, not in % of SSB_max

recr1WithAlleeSorted = exp(y1(indexOfSSbSorted)); %recruitment with Allee

y1WithoutAllee = AnnaModelRecruitmentWithoutAllee(params1);
recr1WithoutAlleeSorted = exp(y1WithoutAllee(indexOfSSbSorted)); %recrutiment without Allee


recr2WithAlleeSorted = exp(y2(indexOfSSbSorted)); %recruitment with Allee

y2WithoutAllee = AnnaModelRecruitmentWithoutAllee(params2);
recr2WithoutAlleeSorted = exp(y2WithoutAllee(indexOfSSbSorted)); %recrutiment without Allee

figure;
hold on;
plot(ssbSorted1, (recr1WithAlleeSorted./ssbSorted), '.-','MarkerSize',10); 
plot(ssbSorted1, (recr1WithoutAlleeSorted./ssbSorted));
xlabel('SSB (% of SSB_{max})');
ylabel('R/SSB');
legend('a = 8.24459539, b = 8.58578E-07, allee = 0.2','a = 8.24459539, b = 8.58578E-07, without Allee');


figure;
hold on;
plot(ssbSorted, (recr2WithAlleeSorted./ssbSorted), '.-','MarkerSize',10); 
plot(ssbSorted, (recr2WithoutAlleeSorted./ssbSorted));
xlabel('SSB (% of SSB_{max})');
ylabel('R/SSB');
legend('a = 8.635711156, b = 2.64626E-06, allee = 0.54595332','a = 8.635711156, b = 2.64626E-06 without Allee');

figure;
hold on;
plot(ssb)

%% Spearman Coefficient

spearmanCorrelationDataAndAnnaModel = corr(y, y1, 'Type', 'Spearman');
