clear
clc

load('behavioralData.mat');

% This script takes a person's blockwise sequence, recodes it to 1 and -1
% and calculates the proportion value for that block. This value is
% averaged over runs and over people. We test for statistical differences 
% between task conditions using ANOVA and multiple comparisons.

% A violin plot is created to show the distribution of proportion values (Fig. 3 in the manuscript)

%% Calculate Proportion value

for sbjct = 1:height(behavioralData)
    sequence = behavioralData.response{sbjct}.keys;
    
    
    % recode HEADS = 1, TAILS = -1
    if string(behavioralData.response{sbjct}.headsPosition(1)) == "right"
        
        sequence(sequence == 1) = -1;    % left outer button, tails, -1
        sequence(sequence == 4) = 1;     % right outer button, heads, 1
        
        sequence(sequence == 2) = -1;    % left inner button, tails, -1
        sequence(sequence == 3) = 1;     % right inner button, heads, 1
        
    else
        
        sequence(sequence == 1) = 1;     % left button, heads, 1 (unnecessary step but left for symmetry)
        sequence(sequence == 4) = -1;    % right button, tails, -1
        
        sequence(sequence == 2) = 1;     % left inner button, heads, 1
        sequence(sequence == 3) = -1;    % right inner button, tails, -1
        
    end
    
    
    for run = 1:6
        
        for block = 1:3
            
            % Extract relevant numbers
            blockwise_array = sequence(behavioralData.response{sbjct}.runNumber == run ...
                & behavioralData.response{sbjct}.blockNumber == block);
            
            % remove NaN values
            blockwise_array(isnan(blockwise_array))= [];
            
            
            % Calculate proportion values
            nHeads = numel(find(blockwise_array == 1));
            nTails = numel(find(blockwise_array == -1));
            prop_perBlock(block) = max(nHeads, nTails)/sum([nHeads,nTails]);
            
            
            
            % calculate mean within block
            mean_perBlock(block) = mean(prop_perBlock);
            
        end
        
        mean_perRun(run) = mean(mean_perBlock);
    end
    mean_perRun_perPerson(sbjct,:) = mean_perRun;
    mean_perPerson(sbjct) = mean(mean_perRun);
   
    
end


%% Average proportion values per condition
mean_proportionValue_perCondition = splitapply(@mean, mean_perPerson', behavioralData.condition);
std_proportionValue_perCondition = splitapply(@std, mean_perPerson', behavioralData.condition);


%% Significance tests -> Is there a significant difference across conditions?
[P,ANOVATAB,STATS]= anova1(mean_perPerson, behavioralData.condition)
[c,m,h,gnames] = multcompare(STATS)

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
    vplot(1,i).MedianColor = "white"
    vplot(1,i).BoxWidth = 0.02;
end
xticklabels(conditionNames)
ylabel("proportion value")
xlabel("condition")
