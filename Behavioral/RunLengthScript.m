clear
clc

load('behavioralData.mat');

% This script calculates the run lengths within every block of every run of every person. 
% The mean across the blocks and across the runs is calculated.
% We test for statistical differences between task conditions using ANOVA .
% A violin plot is created to show the distribution of the average run lengths (Fig. 3 in the manuscript)

%% Calculate average run lengths

for sbjct = 1:height(behavioralData)
    
    % Get subject's sequence
    sequence = behavioralData.response{sbjct}.keys;
    
    for run = 1:6
        for block = 1:3
            
            blockwise_sequence = sequence(behavioralData.response{sbjct}.runNumber == run ...
                & behavioralData.response{sbjct}.blockNumber == block);
            
            % remove NaN values
            blockwise_sequence (isnan(blockwise_sequence ))= [];
            
            % Get length of runs in blockwise_sequence:
            runLengths = diff(find([1,diff(blockwise_sequence'),1]));
            
            % Exception for those people who pressed only one button
            % (irrelevant for current dataset)
            if length(runLengths) > 1
                runLengths(end) = []; % Remove last element because we don't know how long this run could have been
            else
                runLengths = runLengths; % unnecessary step but better readability
            end
            
            % Calculate mean of run lengths in this block
            meanRunLength_perBlock(block) = mean(runLengths);
      
        end
        
        % Collect all mean  values from 6 runs à 3 blocks in one:
        mean_RunLength_perRun(run) = mean(meanRunLength_perBlock);
        
    end
    
    % Collect for every subject:
    mean_RunLength_perRun_perPerson(sbjct,:) = mean_RunLength_perRun;
    mean_RunLength_perPerson(sbjct) = mean(mean_RunLength_perRun);
    
end

mean_RunLength_perPerson = mean_RunLength_perPerson'

%% Average run length values per condition
mean_runLengths_perCondition = splitapply(@mean, mean_RunLength_perPerson, behavioralData.condition);
sd_runLengths_perCondition = splitapply(@std, mean_RunLength_perPerson, behavioralData.condition);

%% Significance test
[P,ANOVATAB,STATS]= anova1(mean_RunLength_perPerson', behavioralData.condition)

%% Fig. 3 Violin plot

for i=1:3
    c = behavioralData.conditionName(behavioralData.condition == i);
    conditionNames(i) = c(1);
end

figure
vplot = violinplot(mean_RunLength_perPerson, behavioralData.condition);
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
xticklabels(conditionNames)
ylabel("Run Length")
