clear; clc;
global recruitment ssb maxSSB lowerSSB higherSSB

year = 2016 %possible values: 2011, 2010
model = 'Anna' %possible values here: 'Anna'
lowerSSB = 0
higherSSB = 0
%% Import Data
if (year == 2010)
    recruitment = xlsread('ourData2010&2011','D3:D49');
    ssb = xlsread('ourData2010&2011','B4:B50');
    maxSSB = max(ssb);
elseif (year == 2011)
    recruitment = xlsread('ourData2010&2011','H3:H50');
    ssb = xlsread('ourData2010&2011','F4:F51');
    maxSSB = max(ssb);
elseif (year == 2016)
    recruitment = xlsread('ourData2010&2011&2016','K3:K55'); % tonnes and thousands
    ssb = xlsread('ourData2010&2011&2016','J4:J56');
    maxSSB = max(ssb);
end
tYears = 1963:year-1;

%% Params for Ricker with Allee

params1 = [1.424005388476614;  -0.000000427484618;   0.265540391475861];
params2 = [1.381793441094838;  -0.000000668699068;   0.098020379944293];
params3 = [1.389953355086599;  -0.000000623498351;   0.199999999999241];

J1 = AnnaModel(params1)
J2 = AnnaModel(params2)
J3 = AnnaModel(params3)

%% Plot graphs Rec(t)

y = log(recruitment);
y1 = AnnaModelRecruitment(params1);
y2 = AnnaModelRecruitment(params2);
y3 = AnnaModelRecruitment(params3);

figure;
hold on;
plot(tYears, y, '*-');
plot(tYears, y1);
legend('data/recruitment', 'a = 1.4240, b = -4.2748E-07, allee = 0.2655');
xlabel('Time(Years)')
ylabel('ln (recruitment)')
title('Comparison ICES 2016 & Ricker With Allee Effect Function');

figure;
hold on;
plot(tYears, y, '*-');
plot(tYears, y2);
legend('data/recruitment', 'a = 1.3817, b = -6.6868E-07, allee = 0.098');
xlabel('Time(Years)')
ylabel('ln (recruitment)')
title('Comparison ICES 2016 & Ricker With Allee Effect Function');

figure;
hold on;
plot(tYears, y, '*-');
plot(tYears, y3);
legend('data/recruitment', 'a = 1.3899, b = -6.2394E-07, allee = 0.1999');
xlabel('Time(Years)')
ylabel('ln (recruitment)')
title('Comparison ICES 2016 & Ricker With Allee Effect Function');

figure;
hold on;
plot(tYears, y, '*-');
plot(tYears, y1, tYears, y2,tYears, y3);
legend('data/recruitment','a = 1.4240, b = -4.2748E-07, allee = 0.2655', 'a = 1.3817, b = -6.6868E-07, allee = 0.098', 'a = 1.3899, b = -6.2394E-07, allee = 0.1999');
xlabel('Time(Years)')
ylabel('ln (recruitment)')
title('Comparison ICES 2016 & Ricker With Allee Effect Function');

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


% recr2WithAlleeSorted = exp(y2(indexOfSSbSorted)); %recruitment with Allee

% y2WithoutAllee = AnnaModelRecruitmentWithoutAllee(params2);
% recr2WithoutAlleeSorted = exp(y2WithoutAllee(indexOfSSbSorted)); %recrutiment without Allee

figure;
hold on;
plot(ssbSorted1, (recr1WithAlleeSorted./ssbSorted), '.-','MarkerSize',10); 
plot(ssbSorted1, (recr1WithoutAlleeSorted./ssbSorted));
xlabel('SSB (% of SSB_{max})');
ylabel('R/SSB');
legend('a = 1.4240, b = -4.2748E-07, allee = 0.2655','a = 1.4240, b = -4.2748E-07, without Allee');

%{
figure;
hold on;
plot(ssbSorted, (recr2WithAlleeSorted./ssbSorted), '.-','MarkerSize',10); 
plot(ssbSorted, (recr2WithoutAlleeSorted./ssbSorted));
xlabel('SSB (% of SSB_{max})');
ylabel('R/SSB');
legend('a = 8.635711156, b = 2.64626E-06, allee = 0.54595332','a = 8.635711156, b = 2.64626E-06 without Allee');
%}
%{
figure;
hold on;
plot(ssb)
%}
%% Spearman Coefficient

[spearmanCorrelationDataAndAnnaModel1, p1] = corr(y, y1, 'Type', 'Spearman')
[spearmanCorrelationDataAndAnnaModel2, p2] = corr(y, y2, 'Type', 'Spearman')
[spearmanCorrelationDataAndAnnaModel3, p3] = corr(y, y3, 'Type', 'Spearman')