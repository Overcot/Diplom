function [xData, x0Data, ssbData, ssbMax, fishMortalityData, gammaData, muData, a, b, allee] = importFromExcel(excelFileName, ourStartingYear, baseStartingYear, T)
    if (strcmp(excelFileName, 'AllDataICES2016'))

        ssbDataRange = strcat('B',int2str(ourStartingYear - baseStartingYear + 2),':','B',int2str(ourStartingYear - baseStartingYear + 2 + T));
        ssbData = xlsread(excelFileName, 1, ssbDataRange);

        xDataRange = strcat('B',int2str(ourStartingYear - baseStartingYear + 2),':','G',int2str(ourStartingYear - baseStartingYear + 2 + T));
        xData = xlsread(excelFileName, 2, xDataRange);
        xData = transpose(xData);
    
        fishMortalityDataRange = strcat('B',int2str(ourStartingYear - baseStartingYear + 2),':','G',int2str(ourStartingYear - baseStartingYear + 2 + T));
        [~, ~, fishMortalityData] = xlsread(excelFileName, 3, fishMortalityDataRange);
        fishMortalityData = str2double(transpose(fishMortalityData));
        
        muDataRange = strcat('A2:F2');
        [~, ~, muData] = xlsread(excelFileName, 4, muDataRange);
        muData = str2double(muData);

        gammaDataRange = strcat('A2:F2');
        [~, ~, gammaData] = xlsread(excelFileName, 5, gammaDataRange);
        gammaData = str2double(gammaData)
        
        x0Data = xData(:,1);
        ssbMax = max(ssbData)
    
        %anna model
        %2016
        %TODO - import this from excel too
        a = 1.389953488259984  
        b = -0.000000623494768   
        allee = 0.199999999966238
        %{
        a = 8.297699765062436 
        b = -0.000000623514772   
        allee = 0.199999967695232
        %}
        %{
        a = 8.28996149679616     
        b = -6.58775962051702e-07
        allee = 0.150063947533343
        %}
        %{
        %2011
        a = 8.244595391;
        b = 8.58578E-07;
        allee = 0.2;
        %}
    else
        error('Not supported import file')
    end

end 