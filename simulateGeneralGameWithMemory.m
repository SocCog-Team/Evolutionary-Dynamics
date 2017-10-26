function [strategyType, strategyStrength] = simulateGeneralGameWithMemory(cfg, nRun, filename)  
  U_shapedDistrCDF = initUshaped();

  
  [pathstr,name,ext] = fileparts(filename) ;  
  if (~isempty(pathstr))
    pathstr = [pathstr '\'];
  end  
  secondaryFilename = [pathstr name '_x' ext];

  
  [~, nStates] = size(cfg.reward);

  
  if ( isfield( cfg, 'propToSeeChoice' ) && isnumeric( cfg.propToSeeChoice ) )
    probabilityToSeeChoice = cfg.propToSeeChoice;
  else
    probabilityToSeeChoice = -1;
  end
  
  if (exist(filename, 'file') && exist(secondaryFilename, 'file'))
    iRun = 0;
    load(filename, 'averageBlindCoopLevel', 'averageOpenCoopLevel', ...
                   'iRun', 'nSurvive', 'meanFit', 'nCooperativeRuns',...
                   'pList', 'qList', 'livingTimeList');  
    load(secondaryFilename, 'xTotalList');                 
    nStartingRun = iRun + 1;             
  else
    averageBlindCoopLevel = zeros(nStates, cfg.nTimeStep, 'single');
    averageOpenCoopLevel = zeros(2*nStates, cfg.nTimeStep, 'single');
    nSurvive = zeros(1, cfg.nTimeStep, 'uint16');
    meanFit  = zeros(1, cfg.nTimeStep);
    nCooperativeRuns = 0;
    pList = [];
    qList = [];
    livingTimeList = [];
    xTotalList = []; 
    nStartingRun = 1;
  end  
  
  for iRun = nStartingRun:nRun   
    nStrategy = 5;
    payoff = zeros(nStrategy);  
    x = (1/nStrategy)*ones(nStrategy, 1);
    livingTime = zeros(nStrategy, 1);
    xTotal = zeros(nStrategy, 1);

    % r - probability to see the choice of the partner
    % p - probabilities to cooperate if previous state was CC, CD, DC and DD
    % q - probabilities to cooperate if previous state was CC, CD, DC and DD and partner cooperates
    p = 0.5*ones(nStates, nStrategy);
    q = 0.5*ones(2*nStates, nStrategy);
    % the probability of seeing has no bias towards edges but is bounded by 0.5         
    if ((probabilityToSeeChoice >= 0) && (probabilityToSeeChoice <= 0.5))
      r = probabilityToSeeChoice*ones(1, nStrategy);
    else
      r = rand(1, nStrategy)/2; 
    end
    
    %init fitness for the first strategy
    payoff = computeExpectedPayoff(payoff, r, p, q, cfg.reward, 1);
    for i = 2:nStrategy      
      p(:, i) = randUshaped(nStates, U_shapedDistrCDF);
      q(:, i) = randUshaped(2*nStates, U_shapedDistrCDF); 
      payoff = computeExpectedPayoff(payoff, r, p, q, cfg.reward, i);
    end

    for t = 1:cfg.nTimeStep
      nSurvive(t) = nSurvive(t) + nnz(x);
      averBlindCoop = p*x;
      averOpenCoop = q*x;
      averageBlindCoopLevel(:, t) = averageBlindCoopLevel(:, t) + averBlindCoop;
      averageOpenCoopLevel(:, t) = averageOpenCoopLevel(:, t) + averOpenCoop;   
      
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
      if (nnz(extinctIndices) > 0)
        %add extinct strategies that lived long enough to the list
        longevityIndices = (extinctIndices & (livingTime > 10*cfg.mutationRate));
        pList = [pList, p(:, longevityIndices)];
        qList = [qList, q(:, longevityIndices)];
        livingTimeList = [livingTimeList; livingTime(longevityIndices)];
        xTotalList = [xTotalList; xTotal(longevityIndices)];
        
        %eliminate rare strategies
        livingTime(extinctIndices) = 0;
        x(extinctIndices) = 0;
        xTotal(extinctIndices) = 0;
        if (x(end-1) == 0) && (x(end) == 0) %if 2 last strategy cells are empty
          %remove the last strategy cell
          x(end) = [];
          livingTime(end) = [];
          xTotal(end) = [];
          r(end) = [];
          p(:, end) = [];
          q(:, end) = [];
          payoff = payoff(1:end-1, 1:end-1);          
          nStrategy = nStrategy - 1;
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
      
        % the probability of seeing has no bias towards edges but is bounded by 0.5 
        if ((probabilityToSeeChoice >= 0) && (probabilityToSeeChoice <= 0.5))
          r(newIndex) = probabilityToSeeChoice;
        else
          r(newIndex) = rand(1)/2; 
        end
        %transform random to U-shaped         
        p(:, newIndex) = randUshaped(nStates, U_shapedDistrCDF);
        q(:, newIndex) = randUshaped(2*nStates, U_shapedDistrCDF);       

        payoff = computeExpectedPayoff(payoff, r, p, q, cfg.reward, newIndex);
        %payoff(x>0, x>0)
      end  
      %re-normalize
      x = x/sum(x); 
    end
    %add remaining strategies to the list
    longevityIndices = (x > 0) & (livingTime > cfg.mutationRate);
    pList = [pList, p(:, longevityIndices)];
    qList = [qList, q(:, longevityIndices)];
    livingTimeList = [livingTimeList; livingTime(longevityIndices)];
    xTotalList = [xTotalList; xTotal(longevityIndices)];
    
    if (totalFit > cfg.coopLowThreshold)
      nCooperativeRuns = nCooperativeRuns + 1;
    end  
    %disp(x);
    %x(x > cfg.initFreq)
    %p(:, x > cfg.initFreq)'
    %q(:, x > cfg.initFreq)'
    save(filename, '-v7.3', 'averageBlindCoopLevel', 'averageOpenCoopLevel', ...
                   'iRun', 'nSurvive', 'meanFit', 'nCooperativeRuns',...
                   'pList', 'qList', 'livingTimeList');  
    save(secondaryFilename, '-v7.3', 'xTotalList');                  
    %disp(q);
  end  
  [strategyType, strategyStrength] = classify_strategy(pList, qList, xTotalList, 11);
  strategyStrength = strategyStrength / (cfg.nTimeStep*nRun/100);
  save(filename, '-v7.3', 'strategyType', 'strategyStrength','-append');            
end
%axis([0 nTimeStep 0 0.01]);


function payoff = computeExpectedPayoff(payoff, r, p, q, reward, currStrategy)
%currStrategy = 1
  [nState, nStrategy] = size(p);

  bothBlindM = zeros(nState);
  firstSeesM = zeros(nState);
  secondSeesM = zeros(nState);
  
  p1Vector = p(:, currStrategy);
  p1Dual = ones(nState, 1) - p1Vector;
  qA1Vector = q(1:nState, currStrategy);
  qA1Dual = ones(nState, 1) - qA1Vector;    
  qB1Vector = q(nState+1:2*nState, currStrategy);
  qB1Dual = ones(nState, 1) - qB1Vector;  
  
  for iStrategy = 1:nStrategy
    %blind actions:
    p2Vector = p(:, iStrategy);
    p2Vector(2) = p(3, iStrategy);
    p2Vector(3) = p(2, iStrategy);
    p2Dual = ones(4, 1) - p2Vector;    
    %if other chooses A
    qA2Vector = q(1:nState, iStrategy);
    qA2Vector(2) = q(3, iStrategy);
    qA2Vector(3) = q(2, iStrategy);
    qA2Dual = ones(nState, 1) - qA2Vector;  
    %if other chooses B   
    qB2Vector = q(nState+1:2*nState, iStrategy);
    qB2Vector(2) = q(7, iStrategy);
    qB2Vector(3) = q(6, iStrategy);
    qB2Dual = ones(nState, 1) - qB2Vector;      

    %bothBlindM - both players do not see actions of the partner;    
    bothBlindM(:, 1) = p1Vector.*p2Vector;
    bothBlindM(:, 2) = p1Vector.*p2Dual;
    bothBlindM(:, 3) = p1Dual.*p2Vector;
    bothBlindM(:, 4) = p1Dual.*p2Dual;
    
    %firstSeesM - first player sees actions of the partner;
    firstSeesM(:, 1) = qA1Vector.*p2Vector;
    firstSeesM(:, 2) = qB1Vector.*p2Dual;
    firstSeesM(:, 3) = qA1Dual.*p2Vector;
    firstSeesM(:, 4) = qB1Dual.*p2Dual;
    
    %secondSeesM - second player sees actions of the partner;
    secondSeesM(:, 1) = p1Vector.*qA2Vector;
    secondSeesM(:, 2) = p1Vector.*qA2Dual;
    secondSeesM(:, 3) = p1Dual.*qB2Vector;
    secondSeesM(:, 4) = p1Dual.*qB2Dual;
    
    M = (1 - r(currStrategy) - r(iStrategy))*bothBlindM + r(currStrategy)*firstSeesM + r(iStrategy)*secondSeesM;
  
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
  