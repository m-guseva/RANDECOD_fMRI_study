function randSeq = createPseudoRandSeq(behavioralData)

% number of trials to produce:
n_trial = height(behavioralData.response{1})
% number of "subjects" (average amount of people per condition)
t = tabulate(behavioralData.condition);
lengthRandSeq = round(mean(t(:,2)));

%Set seed for reproducibility
rng(0,'twister')

% Create pseudorandom sequences
for i = 1:lengthRandSeq
    randSeq{i} = randi([0,1],n_trial,1);
    randSeq{i}(randSeq{i} == 0) = -1; %recode 0 to -1   
end

randSeq = cell2mat(randSeq)'

end