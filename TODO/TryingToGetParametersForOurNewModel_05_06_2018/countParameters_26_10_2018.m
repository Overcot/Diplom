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

%% Try finding parameters using fmincon
% Init data
xdata = ssb;
ydata = log(recruitment);
optionslsqnonlin = optimoptions(@lsqnonlin, 'Algorithm','trust-region-reflective', 'Display', 'iter', 'Tolx',1e-20,'TolFun',1e-20);
A = [];
b = [];
Aeq = [];
beq = [];
for leftABorder=-10000000:0.1:0
    for leftBBorder = leftABorder:0.1:0
        for rightABorder=0:0.1:100000000
            for rightBBorder=0:0.1:10000000
                for leftAllee=0:0.01:1
                    for rightAllee = leftAllee:0.01:1
                        if (model == 'Anna')
                            init = [8,0,0.14];
                            lowerSSB = 0;
                            higherSSB = 0;
                            [x, resnorm] = fmincon(@AnnaModel, init, A,b,Aeq,beq,[leftABorder,leftBBorder,leftAllee], [rightABorder,rightBBorder,rightAllee]);

                        end
                        fileName = 'possibleVariants.xlsx';
                        fileExist = exist(fileName,'file'); 
                        if fileExist==0
                            header = {'A','B','Allee','resnorm', 'leftA','leftB','leftAllee', 'rightA', 'rightB', 'rightAllee'};
                            xlswrite(fileName,header);
                        end
                        [~,~,input] = xlsread(fileName); % Read in your xls file to a cell array (input)
                        new_data = {x(1), x(2), x(3), resnorm, leftABorder, leftBBorder, leftAllee, rightABorder, rightBBorder, rightAllee}; % This is a cell array of the new line you want to add
                        output = cat(1,input,new_data); % Concatinate your new data to the bottom of input
                        xlswrite(fileName,output); % Write to the new excel file. 
                    end
                end
            end
        end
    end
end

