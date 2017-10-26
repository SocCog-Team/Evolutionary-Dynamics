% RT stands for reaction time
function [strategyType, strategyStrength] = simulateIPDwithMemory(cfg, nRun, filename)  
  U_shapedDistrCDF = initUshaped();
  minRT = 2.0;
  maxRT = 3.0;
  stdRT = 0.2;
    
  [pathstr,name,ext] = fileparts(filename) ;  
  if (~isempty(pathstr))
    pathstr = [pathstr '\'];
  end  
  secondaryFilename = [pathstr name '_x' ext];
  
  [~, nStates] = size(cfg.reward);
  
  if ( isfield( cfg, 'probToSeeChoice' ) && isnumeric( cfg.probToSeeChoice ) )
    dltRT = 0.0;
  else
    cfg.probToSeeChoice = -1;
    dltRT = cfg.dltReactionTime;
  end
  
  if (exist(filename, 'file') && exist(secondaryFilename, 'file'))
    iRun = 0;
    load(filename, 'averageBlindCoopLevel', 'averageOpenCoopLevel', 'averageReactionTime', ...
                   'iRun', 'nSurvive', 'meanFit', 'nCooperativeRuns',...
                   'pList', 'qList', 'livingTimeList', 'meanRTList');               
    load(secondaryFilename, 'xTotalList');                 
    nStartingRun = iRun + 1;             
  else
    averageBlindCoopLevel = zeros(nStates, cfg.nTimeStep, 'single');
    averageOpenCoopLevel = zeros(nStates, cfg.nTimeStep, 'single');
    averageReactionTime = zeros(1, cfg.nTimeStep, 'single');
    nSurvive = zeros(1, cfg.nTimeStep, 'uint16');
    meanFit  = zeros(1, cfg.nTimeStep);
    nCooperativeRuns = 0;
    pList = [];
    qList = [];
    livingTimeList = [];
    meanRTList = [];
    xTotalList = []; 
    nStartingRun = 1;
  end  
  
  for iRun = nStartingRun:nRun   
    nStrategy = 5;
    payoff = zeros(nStrategy);  
    x = (1/nStrategy)*ones(nStrategy, 1);
    livingTime = zeros(nStrategy, 1);
    xTotal = zeros(nStrategy, 1);

    % p - probabilities to cooperate if previous state was CC, CD, DC and DD
    % q - probabilities to cooperate if previous state was CC, CD, DC and DD and partner cooperates
    p = 0.5*ones(nStates, nStrategy);
    q = 0.5*ones(nStates, nStrategy);

  	meanRT = minRT + (maxRT - minRT)*rand(1, nStrategy);    
    
    %init fitness for the first strategy
    payoff = computeExpectedPayoff(payoff, p, q, cfg.probToSeeChoice, meanRT, stdRT, dltRT, cfg.reward, 1);
    for i = 2:nStrategy
      p(:, i) = randUshaped(nStates, U_shapedDistrCDF);
      q(:, i) = randUshaped(nStates, U_shapedDistrCDF); 
      payoff = computeExpectedPayoff(payoff, p, q, cfg.probToSeeChoice, meanRT, stdRT, dltRT, cfg.reward, i);
    end

    for t = 1:cfg.nTimeStep
      nSurvive(t) = nSurvive(t) + nnz(x);
      averBlindCoop = p*x;
      averOpenCoop = q*x;
      averRT = dot(meanRT, x);      
      averageBlindCoopLevel(:, t) = averageBlindCoopLevel(:, t) + averBlindCoop;
      averageOpenCoopLevel(:, t) = averageOpenCoopLevel(:, t) + averOpenCoop;   
      averageReactionTime(t) = averageReactionTime(t) + averRT;   
      
      fit = payoff*x;  
      totalFit = sum(x.*fit);
      meanFit(t) = meanFit(t) + totalFit;

      %dx = x.*(fit - totalFit)/totalFit;
      %x = x + dx
      x = (x.*fit)/totalFit;
      x = x/sum(x); %normalize
      
      livingTime = livingTime + 1;
      xTotal = xTotal + x;
      %eliminate all rare strategies
      extinctIndices = ((x > 0) & (x < cfg.minFreq));      
      if ((nnz(extinctIndices) > 0) || (t == cfg.nTimeStep))
        if (nnz(extinctIndices) > 0)
          %add extinct strategies that lived long enough to the list
          longevityIndices = (extinctIndices & (livingTime > 10*cfg.mutationRate));
        else
          %add remaining strategies to the list
          longevityIndices = (x > 0) & (livingTime > 10*cfg.mutationRate);
        end  
        pList = [pList, p(:, longevityIndices)];
        qList = [qList, q(:, longevityIndices)];
        livingTimeList = [livingTimeList; livingTime(longevityIndices)];
        meanRTList = [meanRTList, meanRT(longevityIndices)];
        xTotalList = [xTotalList; xTotal(longevityIndices)];
        
        if (nnz(extinctIndices) > 0)
          %eliminate rare strategies
          livingTime(extinctIndices) = 0;
          x(extinctIndices) = 0;
          xTotal(extinctIndices) = 0;
          if (x(end-1) == 0) && (x(end) == 0) %if 2 last strategy cells are empty
            %remove the last strategy cell
            x(end) = [];
            livingTime(end) = [];
            xTotal(end) = [];
            meanRT(end) = [];
            p(:, end) = [];
            q(:, end) = [];            
            payoff = payoff(1:end-1, 1:end-1);          
            nStrategy = nStrategy - 1;
          end
        end
      end      
  
      %adding new mutant  
      if (rand(1) < 1/cfg.mutationRate)                   
        newIndex = find(x == 0, 1, 'first');
        if (isempty(newIndex))  %add new element for storing strategy
          newIndex = length(x) + 1;
          nStrategy = nStrategy + 1;
        end
        x(newIndex) = cfg.initFreq; 
        livingTime(newIndex) = 0; 
        xTotal(newIndex) = cfg.initFreq;
        
        meanRT(newIndex) = minRT + (maxRT - minRT)*rand(1);
        %transform random to U-shaped        
        p(:, newIndex) = randUshaped(nStates, U_shapedDistrCDF);
        q(:, newIndex) = randUshaped(nStates, U_shapedDistrCDF);        

        payoff = computeExpectedPayoff(payoff, p, q, cfg.probToSeeChoice, meanRT, stdRT, dltRT, cfg.reward, newIndex);
        %payoff(x>0, x>0)
      end  
      %re-normalize
      x = x/sum(x); 
    end    
    if (totalFit > cfg.coopLowThreshold)
      nCooperativeRuns = nCooperativeRuns + 1;
    end  
    %disp(x);
    %x(x > cfg.initFreq)
    %p(:, x > cfg.initFreq)'
    %q(:, x > cfg.initFreq)'
    save(filename, '-v7.3', 'averageBlindCoopLevel', 'averageOpenCoopLevel', 'averageReactionTime', ...
                   'iRun', 'nSurvive', 'meanFit', 'nCooperativeRuns',...
                   'pList', 'qList', 'livingTimeList', 'meanRTList');  
    save(secondaryFilename, '-v7.3', 'xTotalList');                  
    %disp(q);
  end  
  [strategyType, strategyStrength] = classify_strategy(pList, qList, xTotalList, 21);
  strategyStrength = strategyStrength / (cfg.nTimeStep*nRun/100);
  save(filename, '-v7.3', 'strategyType', 'strategyStrength','-append');            
end
%axis([0 nTimeStep 0 0.01]);


function payoff = computeExpectedPayoff(payoff, p, q, probToSeeChoice, meanRT, stdRT, dltRT, reward, currStrategy)
  nStrategy = length(meanRT);

  bothBlindM = zeros(4);
  firstSeesM = zeros(4);
  secondSeesM = zeros(4);
  
  p1Vector = p(:, currStrategy);
  p1Dual = ones(4, 1) - p(:, currStrategy);
  q1Vector = q(:, currStrategy);
  q1Dual = ones(4, 1) - q(:, currStrategy);    
  
  difStdRT = sqrt(2)*stdRT;
  probFirstToSeeChoice = probToSeeChoice; 
  probSecondToSeeChoice = probToSeeChoice;   
 
  for iStrategy = 1:nStrategy %currStrategy is the first player, iStrategy is the second
    if (probToSeeChoice < 0)
      difMeanRT = meanRT(currStrategy) - meanRT(iStrategy); % mean delay of first after second  		
      probFirstToSeeChoice = 1 - normcdf(dltRT, difMeanRT, difStdRT); % if delay > dlt, first sees second;
      probSecondToSeeChoice = normcdf(-dltRT, difMeanRT, difStdRT); % if delay < -dlt, second sees first
    end
    
    p2Vector = p(:, iStrategy);
    p2Vector(2) = p(3, iStrategy);
    p2Vector(3) = p(2, iStrategy);
    p2Dual = ones(4, 1) - p2Vector;
    q2Vector = q(:, iStrategy);
    q2Vector(2) = q(3, iStrategy);
    q2Vector(3) = q(2, iStrategy);
    q2Dual = ones(4, 1) - q2Vector;    

    %bothBlindM - both players do not see actions of the partner;    
    bothBlindM(:, 1) = p1Vector.*p2Vector;
    bothBlindM(:, 2) = p1Vector.*p2Dual;
    bothBlindM(:, 3) = p1Dual.*p2Vector;
    bothBlindM(:, 4) = p1Dual.*p2Dual;
    
    %firstSeesM - first player sees actions of the partner;
    firstSeesM(:, 1) = q1Vector.*p2Vector;
    firstSeesM(:, 2) = 0;
    firstSeesM(:, 3) = q1Dual.*p2Vector;
    firstSeesM(:, 4) = p2Dual;
    
    %secondSeesM - second player sees actions of the partner;
    secondSeesM(:, 1) = p1Vector.*q2Vector;
    secondSeesM(:, 2) = p1Vector.*q2Dual;
    secondSeesM(:, 3) = 0;
    secondSeesM(:, 4) = p1Dual;
    
    M = (1 - probFirstToSeeChoice - probSecondToSeeChoice)*bothBlindM ...
    		+ probFirstToSeeChoice*firstSeesM + probSecondToSeeChoice*secondSeesM;
  
    [statDistr, ~] = eigs( M.', 1 ); %find first left eigenvector
    statDistr = statDistr/sum(statDistr);%normalize it to get stationary distribution 
    payoff(currStrategy, iStrategy) = dot(statDistr, reward(1, :)); 
    payoff(iStrategy, currStrategy) = dot(statDistr, reward(2, :)); 
  end  
end

function U_shapedDistrCDF = initUshaped()
  %generate U-shape distribution
  distrStep = 0.001;
  iDistr = distrStep:distrStep:1-distrStep;   
  U_shapedDistrPDF = 1./(pi*sqrt((1-iDistr).*iDistr));  
  U_shapedDistrCDF = cumsum(U_shapedDistrPDF + [U_shapedDistrPDF(2:end) 0])*distrStep/2;
  %U_shapedDistrCDF = cumsum(U_shapedDistrPDF)*distrStep;
  U_shapedDistrCDF(end) = 1;
end  
  

function randVector = randUshaped(n, U_shapedDistrCDF)
  distrStep = 1/length(U_shapedDistrCDF);
  randValue = rand(1, n);
  randVector = zeros(1, n);
  for i = 1:n
    randVector(i) = (find(U_shapedDistrCDF > randValue(i), 1, 'first') - 1)*distrStep;
    if (randVector(i) < 0.001) 
      randVector(i) = 0.001;      
    end  
    if (randVector(i) > 0.999) 
      randVector(i) = 0.999;      
    end     
  end
end  

%{
  mainIndex = (livingTimeList > 5000000);
  livingTimeList(mainIndex)
  pList(:, mainIndex)
  qList(:, mainIndex)
%}
  