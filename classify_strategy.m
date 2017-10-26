function [sortedStrategyType, sortedStrategyLivingTime] ...
      = classify_strategy(pList, qList, livingTimeList, nLevel)  
    
  levelMark = (0:(nLevel-1))/(nLevel-1); 
  for iLevel = 1:nLevel  
    pList((pList < iLevel/nLevel) & (pList >= (iLevel-1)/nLevel)) = levelMark(iLevel);
    qList((qList < iLevel/nLevel) & (qList >= (iLevel-1)/nLevel)) = levelMark(iLevel);
  end
  strategyList = [pList; qList]';
  
  [strategyType,indexUniqueStrategy,indexStrategy] = unique(strategyList, 'rows');
  nStrategy = length(indexUniqueStrategy); 
  strategyLivingTime = zeros(1, nStrategy);
  for iStrategy = 1:nStrategy
    strategyLivingTime(iStrategy) = sum(livingTimeList(indexStrategy == iStrategy));
  end 
  [sortedStrategyLivingTime, indexSorted] = sort(strategyLivingTime, 'descend');
  sortedStrategyType = strategyType(indexSorted, :)';
end