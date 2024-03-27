clear
clc

load('behavioralData.mat');

% This script extracts the reaction times and calculates the mean RT within 
% each block and averages it over blocks and runs for each subject.
% We test for statistical differences between task conditions using ANOVA and
% perform multiple comparisons.
% A violin plot is created to show the distribution of the average RTs per condition (Fig. 3)



%% Calculate RTs

for sbjct = 1:height(behavioralData)
    array = behavioralData.response{sbjct}.dur_choice;
    
    for run = 1:6
        
        for block = 1:3
            
            % Extract relevant numbers
            blockwise_array = array(behavioralData.response{sbjct}.runNumber == run ...
                & behavioralData.response{sbjct}.blockNumber == block);
            
            % calculate mean within block
            mean_perBlock(block) = mean(blockwise_array, 'omitnan');
            
        end
        
        mean_perRun(run) = mean(mean_perBlock);
    end
    
    mean_perPerson(sbjct) = mean(mean_perRun);
    
end


%% Summary measures for each condition
mean_reactionTime_perCondition = splitapply(@mean, mean_perPerson', behavioralData.condition);
std_reactionTime_perCondition = splitapply(@std, mean_perPerson', behavioralData.condition);

%% Significance test
[p,t,stats] = anova1(mean_perPerson, behavioralData.condition)
%Post hoc test:
[c,m,h,gnames] = multcompare(stats)

%% Fig. 3 Violin plot
for i=1:3
    c = behavioralData.conditionName(behavioralData.condition == i);
    conditionNames(i) = c(1);
end

figure
vplot = violinplot(mean_perPerson, behavioralData.condition);
for i = 1:3
    vplot(1,i).ViolinPlot.FaceColor = "#D3D3D3";
    vplot(1,i).ScatterPlot.MarkerEdgeColor = "#A9A9A9";
    vplot(1,i).ScatterPlot.MarkerFaceColor = "white";
    vplot(1,i).ViolinPlot.LineWidth = 2;
    vplot(1,i).EdgeColor = "#D3D3D3";
    vplot(1,i).BoxPlot.FaceColor = "#3a3a3a";
    vplot(1,i).ScatterPlot.Marker = '.';
    vplot(1,i).ScatterPlot.MarkerFaceAlpha = 0.3;
    vplot(1,i).BoxWidth = 0.02;
end
title(string('Distribution of mean RT per condition'))
xticklabels(conditionNames)
ylabel("seconds")

