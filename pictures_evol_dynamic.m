FontSize = 14;
LineWidth = 2;


%% draw all shares of cooperation (iPD)
shareCooperate = [0.29, 0.55, 0.70, 0.96, 1.00; ...
                  0.25, 0.40, 0.61, 0.75, 0.98; ...
                  0.17, 0.43, 0.67, 0.89, 0.98; ...
                  0.10, 0.30, 0.55, 0.91, 1.00; ...
                  0.09, 0.24, 0.33, 0.55, 0.75; ...
                  0.15, 0.19, 0.18, 0.24, 0.30];
              
figure
set( axes,'fontsize', FontSize, 'FontName', 'Times');
colormap(flipud(autumn(128)));
%colormap(cool(128));
imagesc(shareCooperate);
h = colorbar;
%set( h, 'YDir', 'reverse' );
%hold on;
%for iPlot = 1:6
%  plot(shareCorrect(iPlot, :), lineType{iPlot}, 'linewidth', LineWidth, 'markersize', 9);
%end  
set(gca, 'YDir','normal', 'fontsize', FontSize, 'FontName','Times', ... %'Interpreter', 'latex', ...
         'XTick', 1:5, 'XTickLabel', {'10^4', '10^5', '10^6', '10^7', '10^8'}, ...
         'YTick', 1:6, 'YTickLabel', {'0.0', '0.1', '0.2', '0.3', '0.4', '0.5'});
xlabel( ' Evolution duration [steps] ', 'fontsize', FontSize, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel( ' Probability to see partner`s choice, $\mathrm{p}_\mathrm{see}$', 'fontsize', FontSize, 'FontName', 'Times', 'Interpreter', 'latex');
box off;
%lHandle = legend('p = 0.0', 'p = 0.1',  'p = 0.2',  'p = 0.3',  'p = 0.4',  'p = 0.5', 'location', 'NorthWest');  
%set(lHandle, 'fontsize', FontSize, 'FontName','Times');
set( gcf, 'PaperUnits','centimeters' );
xSize = 24; ySize = 15;
xLeft = 0; yTop = 0;
set( gcf,'PaperPosition', [ xLeft yTop xSize ySize ] );
print ( '-depsc', '-r300', 'shareCooperativeRunsIPD.eps');

%% draw all shares of cooperation (BoS)
shareCooperate = [0.21, 0.38, 0.40, 0.00, 0.30; ...
                  0.04, 0.23, 0.34, 0.00, 0.33; ...
                  0.01, 0.19, 0.37, 0.00, 0.50; ...
                  0.02, 0.19, 0.40, 0.00, 0.55; ...
                  0.05, 0.12, 0.26, 0.00, 0.43; ...
                  0.18, 0.62, 0.99, 0.00, 1.00];

                
figure
set( axes,'fontsize', FontSize, 'FontName', 'Times');
colormap(flipud(autumn(128)));
%colormap(cool(128));
imagesc(shareCooperate);
h = colorbar;
%set( h, 'YDir', 'reverse' );
%hold on;
%for iPlot = 1:6
%  plot(shareCorrect(iPlot, :), lineType{iPlot}, 'linewidth', LineWidth, 'markersize', 9);
%end  
set(gca, 'YDir','normal', 'fontsize', FontSize, 'FontName','Times', ... %'Interpreter', 'latex', ...
         'XTick', 1:5, 'XTickLabel', {'10^4', '10^5', '10^6', '10^7', '10^8'}, ...
         'YTick', 1:6, 'YTickLabel', {'0.0', '0.1', '0.2', '0.3', '0.4', '0.5'});
xlabel( ' Evolution duration [steps] ', 'fontsize', FontSize, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel( ' Probability to see partner`s choice, $\mathrm{p}_\mathrm{see}$', 'fontsize', FontSize, 'FontName', 'Times', 'Interpreter', 'latex');
box off;
%lHandle = legend('p = 0.0', 'p = 0.1',  'p = 0.2',  'p = 0.3',  'p = 0.4',  'p = 0.5', 'location', 'NorthWest');  
%set(lHandle, 'fontsize', FontSize, 'FontName','Times');
set( gcf, 'PaperUnits','centimeters' );
xSize = 24; ySize = 15;
xLeft = 0; yTop = 0;
set( gcf,'PaperPosition', [ xLeft yTop xSize ySize ] );
print ( '-depsc', '-r300', 'shareCooperativeRunsBoS.eps');


%% draw all rewards (iPD)
filenameBase = 'IPD\length8\evolDyn0';

stateName = {'CC', 'CD', 'DC', 'DD'};
nState = 4;
R_VALUE = 3;
P_VALUE = 1;
T_VALUE = 5;
S_VALUE = 0;
evolSettings.reward = [R_VALUE, S_VALUE, T_VALUE, P_VALUE; ...
                       R_VALUE, T_VALUE, S_VALUE, P_VALUE];                   
evolSettings.nTimeStep = 10^8;
outputLength = 10^6;
nTask = 6;
outputRatio = evolSettings.nTimeStep/outputLength;
xValue = outputRatio:outputRatio:evolSettings.nTimeStep;

outputNSurvive = cell(1, nTask);
outputMeanFit = cell(1, nTask);
outputBlindCoopLevel = cell(1, nTask);
outputOpenCoopLevel = cell(1, nTask);

meanP = zeros(6, 4);
meanQ = zeros(6, 4);
stdP = zeros(6, 4);
stdQ = zeros(6, 4);
for iTask = 1:nTask
  filename = [filenameBase num2str(iTask-1) '.mat'];
  load(filename, 'iRun', 'averageBlindCoopLevel', 'averageOpenCoopLevel', ...
                 'nSurvive', 'meanFit', 'pList', 'qList', 'livingTimeList', ...
                 'nCooperativeRuns');
%{
  filenameExt = [filenameBase num2str(iTask) '_x.mat'];               
  load(filename, 'iRun', 'averageBlindCoopLevel', 'averageOpenCoopLevel', ...
                 'nSurvive', 'meanFit', 'pList', 'qList', 'livingTimeList', ...
                 'nCooperativeRuns');               
%}  
  nRun = iRun;              
  outputNSurvive{iTask} = mean(reshape(double(nSurvive)/nRun, outputRatio,[]));
  outputMeanFit{iTask} = mean(reshape(meanFit/nRun, outputRatio,[]));
  averageBlindCoopLevel = averageBlindCoopLevel/nRun;
  averageOpenCoopLevel = averageOpenCoopLevel/nRun;
  outputBlindCoopLevel{iTask} = zeros(nState, outputLength);
  outputOpenCoopLevel{iTask} = zeros(nState, outputLength);
  for iState = 1:nState
    outputBlindCoopLevel{iTask}(iState, :) = mean(reshape(averageBlindCoopLevel(iState, :), outputRatio,[]));
    outputOpenCoopLevel{iTask}(iState, :) = mean(reshape(averageOpenCoopLevel(iState, :), outputRatio,[]));
  end
  
  [meanP(iTask, :), stdP(iTask, :)] = weighted_mean_std(pList, livingTimeList);
  [meanQ(iTask, :), stdQ(iTask, :)] = weighted_mean_std(qList, livingTimeList);
end

%% plot error bars of average strategies (iPD)
stateLabel = {'after CC', 'after CD', 'after DC', 'after DD'};
testLabel = {'0.0', '0.1', '0.2', '0.3', '0.4', '0.5'};

figure('Name', 'Fixations proportions');
set( axes,'fontsize', FontSize, 'FontName','Times');
maxValue = 1.29;

topPlot = subplot(2, 1, 1);
draw_error_bar(meanP, stdP, stateLabel, testLabel, FontSize, 1.49);
xlabel('$\mathrm{p}_\mathrm{see}$', 'fontsize', FontSize, 'FontName','Times', 'Interpreter', 'latex')
box on;
title('cooperation probability without seeing partner`s action', 'fontsize', FontSize, 'FontName','Times', 'Interpreter', 'latex')
bottomPlot = subplot(2, 1, 2);
draw_error_bar(meanQ, stdQ, stateLabel, testLabel, FontSize, 1.49);
xlabel('$\mathrm{p}_\mathrm{see}$', 'fontsize', FontSize, 'FontName','Times', 'Interpreter', 'latex')
box on;
title('cooperation probability when partner is cooperating', 'fontsize', FontSize, 'FontName','Times', 'Interpreter', 'latex')

topSubplotPos = get(topPlot,'Position');
bottomSubplotPos = get(bottomPlot,'Position');
set(bottomPlot,'Position',[bottomSubplotPos(1), topSubplotPos(2)-topSubplotPos(4)-0.1, bottomSubplotPos(3:end)])

set( gcf, 'PaperUnits','centimeters' );
xSize = 28; ySize = 32;
xLeft = 0; yTop = 0;
set( gcf,'PaperPosition', [ xLeft yTop xSize ySize ] );
print ( '-depsc', '-r300', 'barPlotIPD.eps');
  
%% plot rewards (iPD)
lineType = {'-', '-.', '--', ':', '-', '-'};
colorVectorR = [0.0 0.1 0.2 0.4 0.7 1.0];
colorVectorB = [1.0 0.7 0.5 0.2 0.1 0.0];
colorVectorG = [0.0 0.9 0.0 0.6 0.9 0.0];
lineWidth = [1.6, 2.3, 2.3, 2.8, 1.6, 1.6];
lineColor = [1 - colorVectorR; 0.5*colorVectorG; 1 - colorVectorB];
figure

set( axes,'fontsize', FontSize, 'FontName', 'Times');  
hold on
for iTask = 1:nTask
  plot(xValue, movmean(outputMeanFit{iTask}, 5000), lineType{iTask}, 'Color', lineColor(:, iTask), 'linewidth', lineWidth(iTask));
end  
set( gca, 'fontsize', FontSize, 'FontName', 'Times');
%set(gca, 'fontsize', FontSize, 'FontName','Times', 'XTick', 10^7:10^7:10^8, 'XTickLabel', {'10^4', '10^5', '10^6', '10^7', '10^8'});

hold off;
lHandle = legend('$\mathrm{p}_\mathrm{see} = 0.0$', '$\mathrm{p}_\mathrm{see} = 0.1$',  '$\mathrm{p}_\mathrm{see} = 0.2$',  '$\mathrm{p}_\mathrm{see} = 0.3$',  '$\mathrm{p}_\mathrm{see} = 0.4$',  '$\mathrm{p}_\mathrm{see} = 0.5$', 'location', 'SouthEast');  
set(lHandle, 'fontsize', FontSize, 'FontName','Times', 'Interpreter', 'latex');
%axis( [0, max(max(truePPG), max(xEstimPPG)), -0.1, 2.25] );
%set( gca, 'YTick', [], 'fontsize',  FontSize, 'FontName', 'Times');
xlabel( ' Evolution duration [steps] ', 'fontsize', FontSize, 'FontName', 'Times');

set( gcf, 'PaperUnits','centimeters' );
xSize = 24; ySize = 15;
xLeft = 0; yTop = 0;
set( gcf,'PaperPosition', [ xLeft yTop xSize ySize ] );
print ( '-depsc', '-r300', 'meanfitIPD.eps');
  %{
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
%}




%% draw all rewards (BoS)
filenameBase = 'BoS\length8\evolDyn0';

stateName = {'CC', 'CD', 'DC', 'DD'};
nState = 4;
BOTH_OWN_VALUE = 2;
BOTH_OTHER_VALUE = 1;
LEADER_VALUE = 4;
PARTNER_VALUE = 3;
evolSettings.reward = [BOTH_OWN_VALUE, LEADER_VALUE, PARTNER_VALUE, BOTH_OTHER_VALUE; ...
                       BOTH_OWN_VALUE, PARTNER_VALUE, LEADER_VALUE, BOTH_OTHER_VALUE];
evolSettings.nTimeStep = 10^8;
outputLength = 10^6;
nTask = 6;
outputRatio = evolSettings.nTimeStep/outputLength;
xValue = outputRatio:outputRatio:evolSettings.nTimeStep;

outputNSurvive = cell(1, nTask);
outputMeanFit = cell(1, nTask);
outputBlindCoopLevel = cell(1, nTask);
outputOpenCoopLevel = cell(1, nTask);

meanP = zeros(6, 4);
meanQ = zeros(6, 4);
stdP = zeros(6, 4);
stdQ = zeros(6, 4);
for iTask = 1:nTask
  filename = [filenameBase num2str(iTask-1) '.mat'];
  load(filename, 'iRun', 'averageBlindCoopLevel', 'averageOpenCoopLevel', ...
                 'nSurvive', 'meanFit', 'pList', 'qList', 'livingTimeList', ...
                 'nCooperativeRuns');
  nRun = iRun;              
  outputNSurvive{iTask} = mean(reshape(double(nSurvive)/nRun, outputRatio,[]));
  outputMeanFit{iTask} = mean(reshape(meanFit/nRun, outputRatio,[]));
  averageBlindCoopLevel = averageBlindCoopLevel/nRun;
  averageOpenCoopLevel = averageOpenCoopLevel/nRun;
  outputBlindCoopLevel{iTask} = zeros(nState, outputLength);
  outputOpenCoopLevel{iTask} = zeros(nState, outputLength);
  for iState = 1:nState
    outputBlindCoopLevel{iTask}(iState, :) = mean(reshape(averageBlindCoopLevel(iState, :), outputRatio,[]));
    outputOpenCoopLevel{iTask}(iState, :) = mean(reshape(averageOpenCoopLevel(iState, :), outputRatio,[]));
  end
  
  [meanP(iTask, :), stdP(iTask, :)] = weighted_mean_std(pList, livingTimeList);
  [meanQ(iTask, :), stdQ(iTask, :)] = weighted_mean_std(qList, livingTimeList);
end

%% plot error bars of average strategies (BoS)
stateLabel = {'after OWN-OWN', 'after OWN-other', 'after other-OWN', 'after other-other'};
testLabel = {'0.0', '0.1', '0.2', '0.3', '0.4', '0.5'};

figure('Name', 'Fixations proportions');
set( axes,'fontsize', FontSize, 'FontName','Times');
maxValue = 1.29;

topPlot = subplot(3, 1, 1);
draw_error_bar(meanP, stdP, stateLabel, testLabel, FontSize, 1.49);
xlabel('$\mathrm{p}_\mathrm{see}$', 'fontsize', FontSize, 'FontName','Times', 'Interpreter', 'latex')
box on;
title('probability of choosing OWN without seeing partner`s action', 'fontsize', FontSize, 'FontName','Times', 'Interpreter', 'latex')
middlePlot = subplot(3, 1, 2);
draw_error_bar(meanQ(1:4), stdQ(1:4), stateLabel, testLabel, FontSize, 1.49);
xlabel('$\mathrm{p}_\mathrm{see}$', 'fontsize', FontSize, 'FontName','Times', 'Interpreter', 'latex')
box on;
title('probability of choosing OWN when partner chooses OWN', 'fontsize', FontSize, 'FontName','Times', 'Interpreter', 'latex')
bottomPlot = subplot(3, 1, 3);
draw_error_bar(meanQ(5:8), stdQ(5:8), stateLabel, testLabel, FontSize, 1.49);
xlabel('$\mathrm{p}_\mathrm{see}$', 'fontsize', FontSize, 'FontName','Times', 'Interpreter', 'latex')
box on;
title('probability of choosing OWN when partner chooses OWN', 'fontsize', FontSize, 'FontName','Times', 'Interpreter', 'latex')

topSubplotPos = get(topPlot,'Position');
bottomSubplotPos = get(bottomPlot,'Position');
set(bottomPlot,'Position',[bottomSubplotPos(1), topSubplotPos(2)-topSubplotPos(4)-0.1, bottomSubplotPos(3:end)])

set( gcf, 'PaperUnits','centimeters' );
xSize = 28; ySize = 40;
xLeft = 0; yTop = 0;
set( gcf,'PaperPosition', [ xLeft yTop xSize ySize ] );
print ( '-depsc', '-r300', 'barPlotBoS.eps');

%% plot rewards (BoS)
lineType = {'-', '-.', '--', ':', '-', '-'};
colorVectorR = [0.0 0.1 0.2 0.4 0.7 1.0];
colorVectorB = [1.0 0.7 0.5 0.2 0.1 0.0];
colorVectorG = [0.0 0.9 0.0 0.6 0.9 0.0];
lineWidth = [1.6, 2.3, 2.3, 2.8, 1.6, 1.6];
lineColor = [1 - colorVectorR; 0.5*colorVectorG; 1 - colorVectorB];
figure

set( axes,'fontsize', FontSize, 'FontName', 'Times');  
hold on
for iTask = 1:nTask
  plot(xValue, movmean(outputMeanFit{iTask}, 5000), lineType{iTask}, 'Color', lineColor(:, iTask), 'linewidth', lineWidth(iTask));
end  
set( gca, 'fontsize', FontSize, 'FontName', 'Times');
%set(gca, 'fontsize', FontSize, 'FontName','Times', 'XTick', 10^7:10^7:10^8, 'XTickLabel', {'10^4', '10^5', '10^6', '10^7', '10^8'});

hold off;
lHandle = legend('$\mathrm{p}_\mathrm{see} = 0.0$', '$\mathrm{p}_\mathrm{see} = 0.1$',  '$\mathrm{p}_\mathrm{see} = 0.2$',  '$\mathrm{p}_\mathrm{see} = 0.3$',  '$\mathrm{p}_\mathrm{see} = 0.4$',  '$\mathrm{p}_\mathrm{see} = 0.5$', 'location', 'SouthEast');  
set(lHandle, 'fontsize', FontSize, 'FontName','Times', 'Interpreter', 'latex');
%axis( [0, max(max(truePPG), max(xEstimPPG)), -0.1, 2.25] );
%set( gca, 'YTick', [], 'fontsize',  FontSize, 'FontName', 'Times');
xlabel( ' Evolution duration [steps] ', 'fontsize', FontSize, 'FontName', 'Times');

set( gcf, 'PaperUnits','centimeters' );
xSize = 24; ySize = 15;
xLeft = 0; yTop = 0;
set( gcf,'PaperPosition', [ xLeft yTop xSize ySize ] );
print ( '-depsc', '-r300', 'meanfitBoS.eps');
  %{
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
%}



%% tables of strategies
shareWSLS = [56.5 88.9 93.5 82.2 47.2 11.3];
shareALLC = [ 0.0  0.0  0.0  0.0  0.2  2.3];
shareALLD = [ 0.8  4.5  1.8  1.2  5.3  2.7];
shareGTFT = [33.1  0.0  0.1  0.1  0.9  2.9];
shareTFT =  [ 2.0  0.0  0.0  0.0  0.2  0.3];
shareGTFT = shareGTFT + shareTFT;
shareGWLSL = [0.0  0.0  0.2  0.4  2.7  2.0];
shareFbF  = [ 0.0  0.0  0.2  5.1  6.4  3.7];
shareGFbF = [ 0.0  0.0  0.0  1.6  0.5  1.0];
shareGRIM = [ 0.3  0.5  0.1  0.1  1.1  1.4];
X = [0.2 0.4 0.4];
labels = {'Taxes','Expenses','Profit'};
ax1 = subplot(1,2,1);
pie(ax1,X,labels)
title(ax1,'2012');

Y = [0.24 0.46 0.3];
ax2 = subplot(1,2,2);
pie(ax2,Y,labels)
title(ax2,'2013');
