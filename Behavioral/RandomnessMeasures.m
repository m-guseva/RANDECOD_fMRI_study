clear
clc

load('behavioralData.mat');

% This script takes each person's run and calculates the optimal Markov
% Order and the conditional entropy value (at the 3rd order). 

% The values are averaged over runs and an average across conditions is
% calculated.

% We test for statistical differences between task conditions using the Kruskal-Wallis test.
% A violin plot and histogram are created to show the distribution of the
% optimal Markov orders and conditional entropy values.

%% Randomness measure calculated over individual runs
markovOrder_perRun = [];
conditionalEntropy_perRun = [];

for sbjct = 1:height(behavioralData)
    
    for run = 1:6
        run_index = behavioralData.response{sbjct}.runNumber;
        sequence = behavioralData.response{sbjct}.keys(run_index == run);
        
        % remove NaN values
        sequence(isnan(sequence))= [];        
        
        % recode HEADS = 1, TAILS = -1
            if string(behavioralData.response{sbjct}.headsPosition(1)) == "right"
                
                sequence(sequence == 1) = -1;    % left outer button, tails, -1
                sequence(sequence == 4) = 1;     % right outer button, heads, 1
                
                sequence(sequence == 2) = -1;    % left inner button, tails, -1
                sequence(sequence == 3) = 1;     % right inner button, heads, 1
                
            else
                
                sequence(sequence == 1) = 1;     % left button, heads, 1 (unnecessary step but for symmetry)
                sequence(sequence == 4) = -1;    % right button, tails, -1
                
                sequence(sequence == 2) = 1;     % left inner button, heads, 1
                sequence(sequence == 3) = -1;    % right inner button, tails, -1
                
            end
        
        
        % Calculate randomness measures
        [kEst(run), results] = markovOrderBIC(sequence);
        %Extract conditional entropy values (column "CE" in results cell)
        conditionalEntropy = results.CE;
        %Extract value in third position from conditionalEntropy array (=conditional Entropy of third order)
        cEntropyThirdOrder(run) = conditionalEntropy(3);
        
    end
    
    markovOrder_perRun = [markovOrder_perRun; kEst];
    conditionalEntropy_perRun = [conditionalEntropy_perRun; cEntropyThirdOrder];
end

%% Summary measures 
% Average over runs:
median_markovOrder_overRuns = median(markovOrder_perRun,2);
mean_conditionalEntropy_overRuns = mean(conditionalEntropy_perRun,2);

% Average over runs over condition:
median_markovOrder_perCondition = splitapply(@median, median_markovOrder_overRuns, behavioralData.condition);
iqr_markovOrder_perCondition = splitapply(@iqr, median_markovOrder_overRuns, behavioralData.condition);

median_conditionalEntropy_perCondition = splitapply(@median, mean_conditionalEntropy_overRuns, behavioralData.condition);
iqr_conditionalEntropy_perCondition = splitapply(@iqr, mean_conditionalEntropy_overRuns, behavioralData.condition);


%% Significance tests
[p,tbl,stats] = kruskalwallis(median_markovOrder_overRuns, behavioralData.condition)
[p,tbl,stats] = kruskalwallis(mean_conditionalEntropy_overRuns, behavioralData.condition)


%% Plots
for i=1:3
    c = behavioralData.conditionName(behavioralData.condition == i);
    conditionNames(i) = c(1);
end


figure
t = tiledlayout('flow')
for i = 1:3
nexttile
histogram(median_markovOrder_overRuns(behavioralData.condition == i),5, 'FaceColor', '#D3D3D3')
xlabel(conditionNames(i))
ylim([0,20])
xlim([0,2])
if i == 1
    ylabel("Markov Order")
end
end

nexttile([2,3])
vplot = violinplot(mean_conditionalEntropy_overRuns, behavioralData.condition)
for i = 1:3
    vplot(1,i).ViolinPlot.FaceColor = "#D3D3D3";
    vplot(1,i).ScatterPlot.MarkerEdgeColor = "#A9A9A9";
    vplot(1,i).ScatterPlot.MarkerFaceColor = "white";
    vplot(1,i).ViolinPlot.LineWidth = 2;
    vplot(1,i).EdgeColor = "#D3D3D3";
    vplot(1,i).BoxPlot.FaceColor = "#3a3a3a";
    vplot(1,i).ScatterPlot.Marker = '.';
    vplot(1,i).ScatterPlot.MarkerFaceAlpha = 0.3;
    vplot(1,i).MedianColor = "white"
    vplot(1,i).BoxWidth = 0.02;     
end
ylabel("conditional Entropy")
ylim([0,1])
xticklabels(conditionNames)
