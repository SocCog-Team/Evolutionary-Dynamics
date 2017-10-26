  FontSize = 18;
  LineWidth = 1.2;

%1. strategies with memory, various fixed p

%{
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
%}

% iterated BoS
stateName = {'EachOwn', 'Profit1', 'Profit2', 'EachOther'};
nState = 4;
BOTH_OWN_VALUE = 2;
BOTH_OTHER_VALUE = 1;
LEADER_VALUE = 4;
PARTNER_VALUE = 3;
evolSettings.reward = [BOTH_OWN_VALUE, LEADER_VALUE, PARTNER_VALUE, BOTH_OTHER_VALUE; ...
                       BOTH_OWN_VALUE, PARTNER_VALUE, LEADER_VALUE, BOTH_OTHER_VALUE];
evolSettings.coopLowThreshold = 0.95*(LEADER_VALUE + PARTNER_VALUE)/2;

evolSettings.nTimeStep = 10^7;
evolSettings.mutationRate = 100;
evolSettings.minFreq = 0.001; %minimal frequency causing elimination
evolSettings.initFreq = 1.1*evolSettings.minFreq;

outputLength = 5*10^3;
outputRatio = evolSettings.nTimeStep/outputLength;

nRun = 100;
shareCooperativeRuns = zeros(1, 6);

for iTask = 0:5
  filename = ['evolDyn0' num2str(iTask) '.mat'];
  filenameX = ['evolDyn0' num2str(iTask) '_x.mat'];
  evolSettings.propToSeeChoice = iTask/10;

  %output only overaged across all trials
  %to maky plots for a single Run, please call the function for nRun = 1;
  tic
  simulateGeneralGameWithMemory(evolSettings, nRun, filename);
  elapsedTime = toc

  load(filename, 'iRun', 'averageBlindCoopLevel', 'averageOpenCoopLevel', ...
                 'nSurvive', 'meanFit', 'pList', 'qList', 'livingTimeList', ...
                 'nCooperativeRuns');
  nRun = iRun;
  load(filenameX, 'xTotalList');
  if (iTask == 0)
    [w, h] = size(qList);
    qList = zeros(w, h);
  end  
  [strategyType, strategyStrength] = classify_strategy(pList, qList, xTotalList, 21);
  strategyStrength = strategyStrength / (evolSettings.nTimeStep*nRun/100);
  
  nSurvive = double(nSurvive)/nRun;
  meanFit = meanFit/nRun;
  averageBlindCoopLevel = averageBlindCoopLevel/nRun;
  averageOpenCoopLevel = averageOpenCoopLevel/nRun;
  shareCooperativeRuns(iTask + 1) = nCooperativeRuns/nRun;

  outputNSurvive = mean(reshape(nSurvive, outputRatio,[]));
  outputMeanFit = mean(reshape(meanFit, outputRatio,[]));
  outputBlindCoopLevel = zeros(nState, outputLength);
  outputOpenCoopLevel = zeros(2*nState, outputLength);
  for iState = 1:nState
    outputBlindCoopLevel(iState, :) = mean(reshape(averageBlindCoopLevel(iState, :), outputRatio,[]));
    outputOpenCoopLevel(iState, :) = mean(reshape(averageOpenCoopLevel(iState, :), outputRatio,[]));
    outputOpenCoopLevel(iState + nState, :) = mean(reshape(averageOpenCoopLevel(iState + nState, :), outputRatio,[]));
  end
  
    
  disp(strategyStrength(strategyStrength > 1));
  disp(strategyType(:, strategyStrength > 1));
  
  figure
  set( axes,'fontsize', FontSize, 'FontName', 'Times');  
  subplot(3, 2, 1);    
  plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputNSurvive);
  set( gca, 'fontsize', FontSize, 'FontName', 'Times');
  title('Number of strategies', 'fontsize', FontSize, 'FontName', 'Times');
  box off;
    
  subplot(3, 2, 2);    
  plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputMeanFit);
  set( gca, 'fontsize', FontSize, 'FontName', 'Times');
  title('Average reward', 'fontsize', FontSize, 'FontName', 'Times');
  box off;
  
  for iState = 1:nState
    subplot(3, 2, iState + 2);    
    hold on
    plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputBlindCoopLevel(iState, :), 'b-', 'linewidth', LineWidth);
    plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputOpenCoopLevel(iState, :), 'r--', 'linewidth', LineWidth);
    plot(outputRatio:outputRatio:evolSettings.nTimeStep, outputOpenCoopLevel(iState + nState, :), 'm--', 'linewidth', LineWidth);
    hold off
    set( gca, 'fontsize', FontSize, 'FontName', 'Times');
    title(['Average frequency of choosing OWN ' stateName{iState}], 'fontsize', FontSize, 'FontName', 'Times');
    legendHandle = legend('blindly', 'when partner chooses OWN target', 'when partner chooses other target', 'location', 'NorthWest');  
    set(legendHandle, 'fontsize', FontSize, 'FontName', 'Times', 'Interpreter', 'latex');
    box off;
  end  
  set( gcf, 'PaperUnits','centimeters' );
  xSize = 35; ySize = 50;
  xLeft = 0; yTop = 0;
  set( gcf,'PaperPosition', [ xLeft yTop xSize ySize ] );
  print ( '-depsc', '-r300', ['evolDyn0' num2str(iTask) '.eps']);
end



%{
%% rough analysis
share1101 = zeros(1, 6);
share11X1 = zeros(1, 6);
shareCTT = zeros(1, 6);
shareX01X = zeros(1, 6);
share111Xw = zeros(1, 6);
share111X = zeros(1, 6);

evolSettings.nTimeStep = 10^8;
averageP = zeros(6, 4);
averageQ = zeros(6, 8);
for iTask = 1:6
  filename = ['BoS\length8\evolDyn0' num2str(iTask-1) '.mat'];
  filenameX = ['BoS\length8\evolDyn0' num2str(iTask-1) '_x.mat'];
  load(filename, 'iRun', 'averageBlindCoopLevel', 'averageOpenCoopLevel', ...
                 'nSurvive', 'meanFit', 'pList', 'qList', 'livingTimeList', ...
                 'nCooperativeRuns');
  nRun = iRun;
  load(filenameX, 'xTotalList');
                 
  shareCooperativeRuns(iTask) = nCooperativeRuns/nRun;
  
  averageP(iTask, :) = pList*livingTimeList/sum(livingTimeList);
  averageQ(iTask, :) = qList*livingTimeList/sum(livingTimeList);
 
  [strategyType, strategyStrength] = classify_strategy(pList, qList, xTotalList, 21);
    
  strategyStrength = strategyStrength / (evolSettings.nTimeStep*nRun/100);
  
  nStrategy = length(strategyStrength);
  strategyName = cell(1, nStrategy);
  for iStrategy = 1:nStrategy
    stratStr = '            ';
    for i = 1:3*nState
      if (strategyType(i, iStrategy) <= 0.05)
        stratStr(i) = '0';
      elseif (strategyType(i, iStrategy) <= 1/4)
        if (i <= 4)
          stratStr(i) = 'x';
        else
          stratStr(i) = '0';
        end  
      elseif (strategyType(i, iStrategy) >= 0.95)    
        stratStr(i) = '1';            
      elseif (strategyType(i, iStrategy) >= 3/4)  
        if (i <= 4)
          stratStr(i) = 'y';
        else
          stratStr(i) = '1';
        end          
      else
        stratStr(i) = '5';
      end  
    end   
    strategyName{iStrategy} = stratStr;
  end       
  [finalStrategyList,indexUniqueStrategy,indexStrategy] = unique(strategyName);
  nFinalStrategy = length(indexUniqueStrategy); 
  finalStrategyStrength = zeros(1, nFinalStrategy);  
  
  for iStrategy = 1:nFinalStrategy
    finalStrategyStrength(iStrategy) = sum(strategyStrength(indexStrategy == iStrategy));
    if (strcmp(finalStrategyList{iStrategy}(1:3), '110'))
      share1101(iTask) = share1101(iTask) + finalStrategyStrength(iStrategy);
    elseif (strcmp(finalStrategyList{iStrategy}(1:3), '11x'))
      share11X1(iTask) = share11X1(iTask) + finalStrategyStrength(iStrategy);      
    elseif (strcmp(finalStrategyList{iStrategy}(1:4), '5015'))
      shareCTT(iTask) = shareCTT(iTask) + finalStrategyStrength(iStrategy);
    elseif (strcmp(finalStrategyList{iStrategy}(2:3), '01'))
      shareX01X(iTask) = shareX01X(iTask) + finalStrategyStrength(iStrategy);   
    elseif (strcmp(finalStrategyList{iStrategy}(1:3), '111') && ...
            strcmp(finalStrategyList{iStrategy}(5:7), '000') && ...
            strcmp(finalStrategyList{iStrategy}(9:11), '111'))
      share111X(iTask) = share111X(iTask) + finalStrategyStrength(iStrategy);   
    elseif (strcmp(finalStrategyList{iStrategy}(2:3), '11') && ...
            strcmp(finalStrategyList{iStrategy}(6:7), '00'))
      share111Xw(iTask) = share111Xw(iTask) + finalStrategyStrength(iStrategy);       
    end 
  end   
  [sortedStrategyStrength, indexSorted] = sort(finalStrategyStrength, 'descend');
  sortedStrategyList = finalStrategyList(indexSorted);
  
  disp('***********************************');  
  disp(iTask);    
  disp(sortedStrategyStrength(sortedStrategyStrength > 0.25)');
  disp(sortedStrategyList(sortedStrategyStrength > 0.25)');  
end
share1101
share11X1
shareCTT
shareX01X
share111X
share111Xw
%}


 