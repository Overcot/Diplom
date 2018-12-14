global recruitment ssb maxSSB lowerSSB higherSSB

year = 2016 %possible values: 2011, 2010
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
elseif (year == 2016)
    recruitment = xlsread('ourData2010&2011&2016','L3:L55');
    ssb = xlsread('ourData2010&2011&2016','J4:J56');
    maxSSB = max(ssb);
end
%% Try finding parameters using fmincon
% Init data
xdata = ssb;
ydata = log(recruitment);
% options setting
optionslsqnonlin = optimoptions(@lsqnonlin, 'Algorithm','trust-region-reflective', 'Display', 'iter', 'Tolx',1e-20,'TolFun',1e-20);
% u nas net ogranichenii vida Ax <= b & Aeq = B;
A = [];
b = [];
Aeq = [];
beq = [];
% periberaem razlichnie znacheniya dlya levoi i pravoi granici parametrov
% vozmozno(i skorre vsego) oni nevernie

rng default % For reproducibility
gs = GlobalSearch('NumTrialPoints',2000);

problem = createOptimProblem('fmincon','x0',[0,0,0],...
    'objective',@AnnaModel,'lb',[-10,-10,0],'ub',[10,10,1]);
[x,resnorm] = run(gs,problem);

params1 = [8.244595391 8.58578E-07 0.2]
J1 = AnnaModel(params1)

params2 = [8.635711156	2.64626E-06 0.54595332]
J2 = AnnaModel(params2)

params3 = x
J3 = AnnaModel(params3)

problem = createOptimProblem('fmincon','x0',[0,0,0],...
    'objective',@AnnaModel,'lb',[-10,-10,0.1],'ub',[10,10,0.2]);
[x,resnorm] = run(gs,problem);
params4 = x
J4 = AnnaModel(params4)

%{
for leftABorder=-10000000:100000:0
    for leftBBorder = leftABorder:100000:0
        for rightABorder=700000:100000:100000000
            for rightBBorder=0:100000:10000000
                for leftAllee=0:0.1:1
                    % save results in "output" matrix - then in xlsx file
                    output = [];
                    fileName = 'possibleVariants.xlsx';
                    fileExist = exist(fileName,'file'); 
                    if fileExist==0
                        header = {'A','B','Allee','resnorm', 'leftA','leftB','leftAllee', 'rightA', 'rightB', 'rightAllee'};
                        xlswrite(fileName,header);
                    end
                    [~,~,input] = xlsread(fileName); % Read in your xls file to a cell array (input)
                    for rightAllee = leftAllee:0.1:1
                        if (model == 'Anna')
                            init = [leftABorder,leftBBorder,leftAllee];
                            % tried to calculate how many times function goes into one or other side(upper or lower)
                            % more or less then allee effect value
                            lowerSSB = 0; 
                            higherSSB = 0;
                            
                            % there it calculates parameters a, b and allee
                            % threshold with applied borders
          
                            [x, resnorm] = fmincon(@AnnaModel, init, A,b,Aeq,beq,[leftABorder,leftBBorder,leftAllee], [rightABorder,rightBBorder,rightAllee]);
                        end
                        new_data = {x(1), x(2), x(3), resnorm, leftABorder, leftBBorder, leftAllee, rightABorder, rightBBorder, rightAllee}; % This is a cell array of the new line you want to add
                        output = cat(1,output,new_data); % Concatinate your new data to the bottom of input
                    end
                    % save results
                    output = cat(1,input,output); % Concatinate your new data to the bottom of input
                    xlswrite(fileName,output); % Write to the new excel file. 
                end
            end
        end
    end
end
%}
