clear
clc

load('behavioralData.mat');

% This script extracts the responses to the two likert scaled questions that appeared at the end of each of the six runs. 

% We calculate: 
% * average responses across runs per person
% * % of people who had good focus and instruction adherence levels
% * the statistical difference between task conditions using the Kruskal-Wallis Test

% We plot the distribution of answers per run as a stacked bar plot. The 
% responses by conditions are plotted as boxplots. 
% This plot corresponds to Figure 2 in the manuscript.


%% Extract responses on focus and instruction adherence
 
for sbjct = 1:height(behavioralData)
for run = 1:6
focusResponse(sbjct, run) = behavioralData.response{sbjct}.focusResponse(behavioralData.response{1}.runNumber == run ...
    & behavioralData.response{1}.trialNumber == 1 & behavioralData.response{1}.blockNumber== 1 );
instructionsResponse(sbjct, run) = behavioralData.response{sbjct}.instructionsResponse(behavioralData.response{1}.runNumber == run ...
    & behavioralData.response{1}.trialNumber == 1 & behavioralData.response{1}.blockNumber== 1 );
end
end


%% Average response across runs per person

median_focusResponse_overRuns = median(focusResponse,2);
median_instructionsResponse_overRuns = median(instructionsResponse,2);

%% Percent of people who answered good focus and good instructions adherence

tab_focus = tabulate(mode(focusResponse,2));
rel_freq_focus = tab_focus(:,3);

tab_instructions = tabulate(mode(instructionsResponse,2));
rel_freq_instructions = tab_instructions(:,3);

sprintf("%.2f percent of people answered 'very focused' to 'focused'", rel_freq_focus(end)+rel_freq_focus(end-1))
sprintf("%.2f percent of people answered 'very closely' to 'closely'", rel_freq_instructions(end)+rel_freq_instructions(end-1))

sprintf("%.2f percent of people answered 'very focused' ", rel_freq_focus(end-1))
sprintf("%.2f percent of people answered 'very closely' ", rel_freq_instructions(end-1))


%% Fig.2 Distribution of responses per condition

% Calculate the count of each rating category for each rating round
for i = 1:6
    categoryCounts_focusResponse(:, i) = histc(focusResponse(:, i), 1:5);
    categoryCounts_instructionsResponse(:, i) = histc(instructionsResponse(:, i), 1:5);
end

% Calculate the percentage of each rating category for each rating round
categoryPercentages_focusResponse = categoryCounts_focusResponse ./ size(focusResponse, 1);
categoryPercentages_instructionsResponse = categoryCounts_instructionsResponse ./ size(instructionsResponse, 1);

% Plot bar plots and histograms
for i=1:3
    c = behavioralData.conditionName(behavioralData.condition == i);
    conditionNames(i) = c(1);
end

close all
subplot(2,2,1)
ba = bar(flip(categoryPercentages_focusResponse*100)','stacked', 'FaceColor','flat');
colorScheme = linspace(0.9,0, length(ba));

for i = 1:length(ba)
    ba(i).CData = [colorScheme(i), colorScheme(i), colorScheme(i)];
end
ylim([0,100])
xlabel('Run');
ylabel('% of participants answered');
title('Focus');

subplot(2,2,2)
ba = bar(flip(categoryPercentages_instructionsResponse*100)','stacked', 'FaceColor','flat');
for i = 1:length(ba)
    ba(i).CData = [colorScheme(i), colorScheme(i), colorScheme(i)];
end
ylim([0,100])
xlabel('Run');
ylabel('% of participants answered');
title('Instructions adherence');

subplot(2,2,3)
boxplot(median_focusResponse_overRuns, behavioralData.condition, 'Colors', 'k', 'MedianStyle', 'line')
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
xticklabels(conditionNames)
ylim([1,5.5])
yticks([1:5])
ylabel("Rating")

subplot(2,2,4)
boxplot(median_instructionsResponse_overRuns, behavioralData.condition, 'Colors', 'k', 'MedianStyle', 'line')
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
xticklabels(conditionNames)
ylim([1,5.5])
yticks([1:5])


%% Significance tests -> Is there a significant difference across conditions?

[p,tbl,stats] = kruskalwallis(median_focusResponse_overRuns, behavioralData.condition)

[p,tbl,stats] = kruskalwallis(median_instructionsResponse_overRuns, behavioralData.condition)
