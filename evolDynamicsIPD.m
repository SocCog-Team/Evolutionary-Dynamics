%% init game
FontSize = 18;
LineWidth = 1.2;

% iterated Prisoner's dilemma
stateName = {'CC', 'CD', 'DC', 'DD'};
nState = 4;
R_VALUE = 3;
P_VALUE = 1;
T_VALUE = 5;
S_VALUE = 0;
evolSettings.reward = [R_VALUE, S_VALUE, T_VALUE, P_VALUE; ...
                       R_VALUE, T_VALUE, S_VALUE, P_VALUE];
evolSettings.coopLowThreshold = 0.9*R_VALUE;                     

evolSettings.nTimeStep = 10^8;
evolSettings.mutationRate = 100;
evolSettings.minFreq = 0.001; %minimal frequency causing elimination
evolSettings.initFreq = 1.1*evolSettings.minFreq;

outputLength = 10^5;
outputRatio = evolSettings.nTimeStep/outputLength;

nRun = 40;
shareCooperativeRuns = zeros(1, 6);
shareGTFT = zeros(1, 6);
shareWSLS = zeros(1, 6);
shareALLC = zeros(1, 6);
shareALLD = zeros(1, 6);
shareTFT = zeros(1, 6);
shareFbF = zeros(1, 6);
shareGFbF = zeros(1, 6);
shareGWLSL = zeros(1, 6);
shareGRIM = zeros(1, 6);


shareTotalGTFT = zeros(1, 6);
shareTotalWSLS = zeros(1, 6);
shareTotalALLC = zeros(1, 6);
shareTotalALLD = zeros(1, 6);
shareTotalTFT = zeros(1, 6);
shareTotalFbF = zeros(1, 6);
shareTotalGFbF = zeros(1, 6);
shareTotalGWLSL = zeros(1, 6);
shareTotalGRIM = zeros(1, 6);

%{
%% run with fixed probability to see opponent choice
for iTask = 0:5
  filename = ['evolDyn0' num2str(iTask) '.mat'];
  evolSettings.probToSeeChoice = iTask/10;

  %output only overaged across all trials
  %to maky plots for a single Run, please call the function for nRun = 1;
  tic
  [strategyType, strategyStrength] = simulateIPDwithMemory(evolSettings, nRun, filename);
  elapsedTime = toc

  load(filename, 'iRun', 'averageBlindCoopLevel', 'averageOpenCoopLevel', ...
                 'nSurvive', 'meanFit', 'pList', 'qList', 'livingTimeList', ...
                 'nCooperativeRuns');
               
  %load(filename, 'strategyType', 'strategyStrength');                  
  nSurvive = double(nSurvive)/nRun;
  meanFit = meanFit/nRun;
  averageBlindCoopLevel = averageBlindCoopLevel/nRun;
  averageOpenCoopLevel = averageOpenCoopLevel/nRun;
  shareCooperativeRuns(iTask + 1) = nCooperativeRuns/nRun;

  outputNSurvive = mean(reshape(nSurvive, outputRatio,[]));
  outputMeanFit = mean(reshape(meanFit, outputRatio,[]));
  outputBlindCoopLevel = zeros(nState, outputLength);
  outputOpenCoopLevel = zeros(nState, outputLength);
  for iState = 1:nState
    outputBlindCoopLevel(iState, :) = mean(reshape(averageBlindCoopLevel(iState, :), outputRatio,[]));
    outputOpenCoopLevel(iState, :) = mean(reshape(averageOpenCoopLevel(iState, :), outputRatio,[]));
  end
  
  
 
  index = find(ismember(strategyType(1:4, :)', [1, 0, 0, 1], 'rows'));
  shareWSLS(iTask + 1) = sum(strategyStrength(index));
  shareGTFT(iTask + 1) = sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0.5, 1, 0.5], 'rows')));
  shareALLD(iTask + 1) = sum(strategyStrength( ismember(strategyType(1:4, :)', [0, 0, 0, 0], 'rows')));
  shareALLC(iTask + 1) = sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 1, 1, 1], 'rows')));
  shareTFT(iTask + 1) = sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0, 1, 0], 'rows')));
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0, 1, 0.5], 'rows'))));
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0.5, 1, 0], 'rows'))));
  disp('GWLSL:');
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0, 0.5, 1], 'rows'))));
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0.5, 0, 1], 'rows'))));
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0.5, 0.5, 1], 'rows'))));

  
  
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [1, 0, 0, 0], 'rows')))));
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [0.5, 0, 0, 0], 'rows')))));
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [1, 0, 1, 1], 'rows')))));
   
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [0, 1, 1, 0.5], 'rows')))));
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [1, 1, 1, 0.5], 'rows')))));
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [1, 1, 1, 0], 'rows')))));
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [0.5, 1, 1, 0], 'rows')))));
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [0, 1, 1, 0], 'rows')))));
    
  disp(strategyStrength(strategyStrength > 1));
  disp(strategyType(:, strategyStrength > 1));
  
  figure
  set( axes,'fontsize', FontSize, 'FontName', 'Times');  
  subplot(nState + 2, 1, 1);    
  plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputNSurvive, 'linewidth', LineWidth);
  set( gca, 'fontsize', FontSize, 'FontName', 'Times');
  box off;
  title('Number of strategies', 'fontsize', FontSize, 'FontName', 'Times');
    
  subplot(nState + 2, 1, 2);    
  plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputMeanFit, 'linewidth', LineWidth);
  set( gca, 'fontsize', FontSize, 'FontName', 'Times');
  box off;
  title('Average reward', 'fontsize', FontSize, 'FontName', 'Times');
  
  for iState = 1:nState
    subplot(nState + 2, 1, iState + 2);    
    hold on
    plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputBlindCoopLevel(iState, :), 'b-', 'linewidth', LineWidth);
    plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputOpenCoopLevel(iState, :), 'r--', 'linewidth', LineWidth);
    hold off
    legend('blindly', 'when partner cooperate', 'location', 'NorthWest');  
    set( gca, 'fontsize', FontSize, 'FontName', 'Times');
    box off;
    title(['Average cooperation rate after ' stateName{iState}], 'fontsize', FontSize, 'FontName', 'Times');
  end  
  set( gcf, 'PaperUnits','centimeters' );
  xSize = 35; ySize = 50;
  xLeft = 0; yTop = 0;
  set( gcf,'PaperPosition', [ xLeft yTop xSize ySize ] );
  print ( '-depsc', '-r300', ['evolDyn0' num2str(iTask) '.eps']);
end
%}


%%{
%% run with fixed probability to see opponent choice
evolSettings.nTimeStep = 10^8;
averageP = zeros(6, 4);
averageQ = zeros(6, 4);
for iTask = 1:6
  filename = ['IPD\length8\evolDyn0' num2str(iTask-1) '.mat'];
  filenameX = ['IPD\length8\evolDyn0' num2str(iTask-1) '_x.mat'];
  load(filename, 'iRun', 'averageBlindCoopLevel', 'averageOpenCoopLevel', ...
                 'nSurvive', 'meanFit', 'pList', 'qList', 'livingTimeList', ...
                 'nCooperativeRuns');
               
  averageP(iTask, :) = pList*livingTimeList/sum(livingTimeList);
  averageQ(iTask, :) = qList*livingTimeList/sum(livingTimeList);
  
  
  if (iTask == 1)
    continue
  end
  load(filenameX, 'xTotalList');  
  [strategyType, strategyStrength] = classify_strategy(pList, qList, xTotalList, 21);
    
  strategyStrength = strategyStrength / (evolSettings.nTimeStep*nRun/100);
  
  nStrategy = length(strategyStrength);
  strategyName = cell(1, nStrategy);
  subStrategyLength = 4;  
  maxThresh = [0.2, 0.25, 0.25];
  sumThresh = [0.4, 0.5, 0.5];
  for iStrategy = 1:nStrategy
    strategyName{iStrategy} = '';
    coarseStrategy = round(strategyType(:, iStrategy));
    stratStr = num2str(coarseStrategy');
    stratStr = stratStr(stratStr ~= ' ');
    stratError = abs(strategyType(:, iStrategy) - coarseStrategy);    
    nSubStrategySection = ceil(length(coarseStrategy)/subStrategyLength);
    for iStrategySection = 1:nSubStrategySection
      index = (1:4) + 4*(iStrategySection - 1);
      subStrategyError = stratError(index);
      if (nnz(subStrategyError) == 0)
        strategyName{iStrategy} = [strategyName{iStrategy}, stratStr(index)];
      elseif ((max(subStrategyError) <= maxThresh(iStrategySection)) && (sum(subStrategyError) <= sumThresh(iStrategySection)))
        strategyName{iStrategy} = [strategyName{iStrategy}, stratStr(index),'w'];
      else
        for i = index
          if (strategyType(i, iStrategy) <= 0.11)
            stratStr(i) = '0';
          elseif (strategyType(i, iStrategy) < 1/3)
            stratStr(i) = 'x';
          elseif (strategyType(i, iStrategy) >= 0.89)    
            stratStr(i) = '1';            
          elseif (strategyType(i, iStrategy) > 2/3)  
            stratStr(i) = 'y';
          else
            stratStr(i) = '5';
          end  
        end  
        strategyName{iStrategy} = [strategyName{iStrategy}, stratStr(index)];
      end    
    end        
    if (strategyName{iStrategy}(subStrategyLength+1) == 'w')
      shift = 6;
    else
      shift = 5;
    end
    [~, isPreciseStrat] = str2num(strategyName{iStrategy}(shift:shift+subStrategyLength-1));
    if (~isPreciseStrat)
      strategyName{iStrategy}(shift:end) = [];
      if (sum(strategyType(end-subStrategyLength+1:end, iStrategy)) > subStrategyLength/2)
        strategyName{iStrategy} = strcat(strategyName{iStrategy}, 'COOP');
      else
        strategyName{iStrategy} = strcat(strategyName{iStrategy}, 'NOCO');
      end  
    end  
  end       
  [finalStrategyList,indexUniqueStrategy,indexStrategy] = unique(strategyName);
  nFinalStrategy = length(indexUniqueStrategy); 
  finalStrategyStrength = zeros(1, nFinalStrategy);  
  
  for iStrategy = 1:nFinalStrategy
    finalStrategyStrength(iStrategy) = sum(strategyStrength(indexStrategy == iStrategy));
    if (strcmp(finalStrategyList{iStrategy}(1:4), '1001'))
      shareWSLS(iTask) = shareWSLS(iTask) + finalStrategyStrength(iStrategy);
    elseif (strcmp(finalStrategyList{iStrategy}(1:4), '1111'))
      shareALLC(iTask) = shareALLC(iTask) + finalStrategyStrength(iStrategy);
    elseif (strcmp(finalStrategyList{iStrategy}(1:4), '0000'))
      shareALLD(iTask) = shareALLD(iTask) + finalStrategyStrength(iStrategy);
    elseif (strcmp(finalStrategyList{iStrategy}(1:4), '1x1x') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), '1x15') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), '1x1y') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), '1y15') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), '1y1y') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), '1y1x') ||...            
            strcmp(finalStrategyList{iStrategy}(1:4), '1515') ||...            
            strcmp(finalStrategyList{iStrategy}(1:4), '151y') ||...            
            strcmp(finalStrategyList{iStrategy}(1:4), '151x'))
      shareGTFT(iTask) = shareGTFT(iTask) + finalStrategyStrength(iStrategy);
    elseif (strcmp(finalStrategyList{iStrategy}(1:4), '1xx1') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), '1x51') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), '1xy1') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), '1y51') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), '1yy1') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), '1yx1') ||...            
            strcmp(finalStrategyList{iStrategy}(1:4), '1551') ||...            
            strcmp(finalStrategyList{iStrategy}(1:4), '15y1') ||...            
            strcmp(finalStrategyList{iStrategy}(1:4), '15x1'))
      shareGWLSL(iTask) = shareGWLSL(iTask) + finalStrategyStrength(iStrategy);      
    elseif (strcmp(finalStrategyList{iStrategy}(1:4), '1010'))
      shareTFT(iTask) = shareTFT(iTask) + finalStrategyStrength(iStrategy);  
    elseif (strcmp(finalStrategyList{iStrategy}(1:4), '1011'))
      shareFbF(iTask) = shareFbF(iTask) + finalStrategyStrength(iStrategy); 
    elseif (strcmp(finalStrategyList{iStrategy}(1:4), '1x11') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), '1511') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), '1y11'))  
      shareGFbF(iTask) = shareGFbF(iTask) + finalStrategyStrength(iStrategy); 
    elseif (strcmp(finalStrategyList{iStrategy}(1:4), '1000') ||...
            strcmp(finalStrategyList{iStrategy}(1:4), 'y000')) 
      shareGRIM(iTask) = shareGRIM(iTask) + finalStrategyStrength(iStrategy); 
    end  
  end 
  [sortedStrategyStrength, indexSorted] = sort(finalStrategyStrength, 'descend');
  sortedStrategyList = finalStrategyList(indexSorted);
  
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0, 1, 0.5], 'rows'))));
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0.5, 1, 0], 'rows'))));
  disp('GWLSL:');
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0, 0.5, 1], 'rows'))));
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0.5, 0, 1], 'rows'))));
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0.5, 0.5, 1], 'rows'))));

  
  
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [1, 0, 0, 0], 'rows')))));
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [0.5, 0, 0, 0], 'rows')))));
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [1, 0, 1, 1], 'rows')))));
   
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [0, 1, 1, 0.5], 'rows')))));
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [1, 1, 1, 0.5], 'rows')))));
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [1, 1, 1, 0], 'rows')))));
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [0.5, 1, 1, 0], 'rows')))));
%  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [0, 1, 1, 0], 'rows')))));
  disp('***********************************');  
  disp(iTask);    
  disp(sortedStrategyStrength(sortedStrategyStrength > 0.25)');
  disp(sortedStrategyList(sortedStrategyStrength > 0.25)');
  
  totalStrategyName = cell(1, nStrategy);
  totalStrategyType = (1 - iTask/10)*strategyType(1:4, :) + (iTask/10)*strategyType(5:8, :)
  for iStrategy = 1:nStrategy   
    coarseStrategy = round(totalStrategyType(:, iStrategy));
    stratStr = num2str(coarseStrategy');
    stratStr = stratStr(stratStr ~= ' ');
    strategyError = abs(totalStrategyType(:, iStrategy) - coarseStrategy);    
    
    if (nnz(stratStr) == 0)
      totalStrategyName{iStrategy} = stratStr;
    elseif ((max(strategyError) <= maxThresh(1)) && (sum(strategyError) <= sumThresh(1)))
      totalStrategyName{iStrategy} = [stratStr,'w'];
    else
      for i = 1:4
        if (strategyType(i, iStrategy) <= 0.11)
          stratStr(i) = '0';
        elseif (strategyType(i, iStrategy) < 1/3)
          stratStr(i) = 'x';
        elseif (strategyType(i, iStrategy) >= 0.89)    
          stratStr(i) = '1';            
        elseif (strategyType(i, iStrategy) > 2/3)  
          stratStr(i) = 'y';
        else
          stratStr(i) = '5';
        end  
      end  
      totalStrategyName{iStrategy} = stratStr;
    end    
  end        
  [finalTotalStrategyList,indexUniqueStrategy,indexStrategy] = unique(totalStrategyName);
  nFinalStrategy = length(indexUniqueStrategy); 
  finalTotalStrategyStrength = zeros(1, nFinalStrategy);  
  
  for iStrategy = 1:nFinalStrategy
    finalTotalStrategyStrength(iStrategy) = sum(strategyStrength(indexStrategy == iStrategy));
    if (strcmp(finalTotalStrategyList{iStrategy}, '1001'))
      shareTotalWSLS(iTask) = shareTotalWSLS(iTask) + finalTotalStrategyStrength(iStrategy);
    elseif (strcmp(finalTotalStrategyList{iStrategy}(1:4), '1111'))
      shareTotalALLC(iTask) = shareTotalALLC(iTask) + finalTotalStrategyStrength(iStrategy);
    elseif (strcmp(finalTotalStrategyList{iStrategy}(1:4), '0000'))
      shareTotalALLD(iTask) = shareTotalALLD(iTask) + finalTotalStrategyStrength(iStrategy);
    elseif (strcmp(finalTotalStrategyList{iStrategy}(1:4), '1x1x') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1x15') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1x1y') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1y15') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1y1y') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1y1x') ||...            
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1515') ||...            
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '151y') ||...            
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '151x'))
      shareTotalGTFT(iTask) = shareTotalGTFT(iTask) + finalTotalStrategyStrength(iStrategy);
    elseif (strcmp(finalTotalStrategyList{iStrategy}(1:4), '1xx1') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1x51') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1xy1') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1y51') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1yy1') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1yx1') ||...            
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1551') ||...            
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '15y1') ||...            
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '15x1'))
      shareTotalGWLSL(iTask) = shareTotalGWLSL(iTask) + finalTotalStrategyStrength(iStrategy);      
    elseif (strcmp(finalTotalStrategyList{iStrategy}(1:4), '1010'))
      shareTotalTFT(iTask) = shareTotalTFT(iTask) + finalTotalStrategyStrength(iStrategy);  
    elseif (strcmp(finalTotalStrategyList{iStrategy}(1:4), '1011'))
      shareTotalFbF(iTask) = shareTotalFbF(iTask) + finalTotalStrategyStrength(iStrategy); 
    elseif (strcmp(finalTotalStrategyList{iStrategy}(1:4), '1x11') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1511') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), '1y11'))  
      shareTotalGFbF(iTask) = shareTotalGFbF(iTask) + finalTotalStrategyStrength(iStrategy); 
    elseif (strcmp(finalTotalStrategyList{iStrategy}(1:4), '1000') ||...
            strcmp(finalTotalStrategyList{iStrategy}(1:4), 'y000')) 
      shareTotalGRIM(iTask) = shareTotalGRIM(iTask) + finalTotalStrategyStrength(iStrategy); 
    end  
  end 
  [sortedTotalStrategyStrength, indexSorted] = sort(finalTotalStrategyStrength, 'descend');
  sortedTotalStrategyList = finalTotalStrategyList(indexSorted); 
end
%%}


%% run with mutable probability to see opponent choice

%Pr to see "conspecifics" action: 0.01 |0.1   |0.2  |0.3   |0.33  |0.4   |0.49
%dltReactionTime                : 0.658|0.3625|0.238|0.1483|0.1218|0.0716|0.007

allDltReactionTime = [0.658, 0.3625, 0.238, 0.1483, 0.0716, 0.007];
for iTask = 1:6
  filename = ['evolDynRT' num2str(iTask) '.mat'];
  evolSettings.dltReactionTime = allDltReactionTime(iTask);

  %output only overaged across all trials
  %to maky plots for a single Run, please call the function for nRun = 1;
  tic
  [strategyType, strategyStrength] = simulateIPDwithMemory(evolSettings, nRun, filename);
  elapsedTime = toc

  load(filename, 'iRun', 'averageBlindCoopLevel', 'averageOpenCoopLevel', 'averageReactionTime',...
                 'nSurvive', 'meanFit', 'pList', 'qList', 'livingTimeList', ...
                 'meanRTList', 'nCooperativeRuns');
               
  %load(filename, 'strategyType', 'strategyStrength');                  
  nSurvive = double(nSurvive)/nRun;
  meanFit = meanFit/nRun;
  averageBlindCoopLevel = averageBlindCoopLevel/nRun;
  averageOpenCoopLevel = averageOpenCoopLevel/nRun;
  averageReactionTime = averageReactionTime/nRun;
  shareCooperativeRuns(iTask) = nCooperativeRuns/nRun;

  outputNSurvive = mean(reshape(nSurvive, outputRatio,[]));
  outputMeanFit = mean(reshape(meanFit, outputRatio,[]));
  outputReactionTime = mean(reshape(averageReactionTime, outputRatio,[]));
  outputBlindCoopLevel = zeros(nState, outputLength);
  outputOpenCoopLevel = zeros(nState, outputLength);
  for iState = 1:nState
    outputBlindCoopLevel(iState, :) = mean(reshape(averageBlindCoopLevel(iState, :), outputRatio,[]));
    outputOpenCoopLevel(iState, :) = mean(reshape(averageOpenCoopLevel(iState, :), outputRatio,[]));
  end
      
 
  index = find(ismember(strategyType(1:4, :)', [1, 0, 0, 1], 'rows'));
  shareWSLS(iTask) = sum(strategyStrength(index));
  shareGTFT(iTask) = sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0.5, 1, 0.5], 'rows')));
  shareALLD(iTask) = sum(strategyStrength( ismember(strategyType(1:4, :)', [0, 0, 0, 0], 'rows')));
  shareALLC(iTask) = sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 1, 1, 1], 'rows')));
  shareTFT(iTask) = sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0, 1, 0], 'rows')));
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0, 1, 0.5], 'rows'))));
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0.5, 1, 0], 'rows'))));
  disp('GWLSL:');
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0, 0.5, 1], 'rows'))));
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0.5, 0, 1], 'rows'))));
  disp(sum(strategyStrength( ismember(strategyType(1:4, :)', [1, 0.5, 0.5, 1], 'rows'))));

  
%{  
  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [1, 0, 0, 0], 'rows')))));
  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [0.5, 0, 0, 0], 'rows')))));
  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [1, 0, 1, 1], 'rows')))));
  
  
  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [0, 1, 1, 0.5], 'rows')))));
  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [1, 1, 1, 0.5], 'rows')))));
  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [1, 1, 1, 0], 'rows')))));
  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [0.5, 1, 1, 0], 'rows')))));
  disp(sum(strategyStrength( find(ismember(strategyType(1:4, :)', [0, 1, 1, 0], 'rows')))));
%}    
  disp(strategyStrength(strategyStrength > 1));
  disp(strategyType(:, strategyStrength > 1));
  
  figure
  set( axes,'fontsize', FontSize, 'FontName', 'Times');  
  subplot(nState, 2, 1);    
  plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputNSurvive, 'linewidth', LineWidth);
  set( gca, 'fontsize', FontSize, 'FontName', 'Times');
  box off;
  title('Number of strategies', 'fontsize', FontSize, 'FontName', 'Times');
    
  subplot(nState, 2, 3);    
  plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputMeanFit, 'linewidth', LineWidth);
  set( gca, 'fontsize', FontSize, 'FontName', 'Times');
  box off;
  title('Average reward', 'fontsize', FontSize, 'FontName', 'Times');

  subplot(nState, 2, 5);    
  plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputReactionTime, 'linewidth', LineWidth);
  set( gca, 'fontsize', FontSize, 'FontName', 'Times');
  box off;
  title('Average reaction time', 'fontsize', FontSize, 'FontName', 'Times');
  
  
  for iState = 1:nState
    subplot(nState, 2, 2*iState);    
    hold on
    plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputBlindCoopLevel(iState, :), 'b-', 'linewidth', LineWidth);
    plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputOpenCoopLevel(iState, :), 'r--', 'linewidth', LineWidth);
    hold off
    legend('blindly', 'when partner cooperate', 'location', 'NorthWest');  
    set( gca, 'fontsize', FontSize, 'FontName', 'Times');
    box off;
    title(['Average cooperation rate after ' stateName{iState}], 'fontsize', FontSize, 'FontName', 'Times');
  end  
  set( gcf, 'PaperUnits','centimeters' );
  xSize = 35; ySize = 50;
  xLeft = 0; yTop = 0;
  set( gcf,'PaperPosition', [ xLeft yTop xSize ySize ] );
  print ( '-depsc', '-r300', ['evolDyn0' num2str(iTask) '.eps']);
end
%{
disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 1, 1, 1], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 0, 0, 1], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 0, 0, 0], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 1, 1, 0], 'rows')))));

disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 0, 1, 0], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 0.5, 0, 1], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 0, 1, 1], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 0.5, 1, 0.5], 'rows')))));

disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 1, 0, 1], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 1, 0.5, 1], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 0, 0.5, 1], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 0.5, 0.5, 1], 'rows')))));

disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 0, 1, 0.5], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 1, 0.5, 1, 1], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType', [1, 0, 0, 1, 0, 0, 0, 0], 'rows')))));




disp(sum(strategyStrength( find(ismember(strategyType(5:8, :)', [0, 0, 0, 0], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType(5:8, :)', [1, 0, 0, 0], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType(5:8, :)', [0.5, 0, 0, 0], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType(5:8, :)', [1, 0, 0, 1], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType(5:8, :)', [1, 1, 1, 1], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType(5:8, :)', [1, 1, 0, 1], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType(5:8, :)', [1, 1, 0.5, 1], 'rows')))));

disp(sum(strategyStrength( find(ismember(strategyType(5:8, :)', [1, 1, 1, 0], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType(5:8, :)', [1, 0, 1, 0], 'rows')))));
disp(sum(strategyStrength( find(ismember(strategyType(5:8, :)', [1, 0, 0.5, 1], 'rows')))));

disp(sum(strategyStrength( find(ismember(strategyType(5:8, :)', [0.5, 0, 0.5, 0 ], 'rows')))));


%}
 