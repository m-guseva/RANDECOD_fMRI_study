function [kEst, results, P1] = markovOrderBIC(s, kMax, segOff, segLen)

% estimation of the order of a Markov process using the BIC criterion
%
% [kEst, results] = markovOrderBIC(s, kMax, segOff, segLen)
%
% s:        sequence of state labels
% kMax:     largest order to consider (may be inf)
% segOff:   offsets of sequence segments, zero-based
% segLen:   length of sequence segments
%
% kEst:     estimated Markov order: the value of k for which BIC is minimal
% results:  table with variables
%           k:          order
%           m2lL:       -2 log likelihood at ML-estimated parameters
%           penalty:    BIC penalty term
%           BIC:        BIC, m2lL + penalty
%           CE:         conditional entropy, in bit
% following Csiszár and Shields (2000), "The consistency of the BIC Markov
% order estimator", The Annals of Statistics 28(6), 1601-1619

% This function was written by Dr. Carsten Allefeld (https://allefeld.github.io)


% transform sequence into column vector
s = s(:);
% length of sequence
n = numel(s);

% if no segment offsets specified, treat sequence as one segment
if ~exist('segOff', 'var')
    segOff = 0;
end

% if no segment lengths specified, derive from offsets and sequence length
if ~exist('segLen', 'var')
    segLen = diff([segOff, n]);
end

% if no maximum order specified or inf, go up to maximum
if ~exist('kMax', 'var') || isinf(kMax)
    kMax = max(segLen) - 1;
end


% translate occurring states to integers
[stateLabels, ~, state] = unique(s);

% number of states
nStates = numel(stateLabels);


% initialize results for requested maximum Markov order
m2lL = zeros(kMax + 1, 1);
CE = zeros(kMax + 1, 1);

% compute m2lL and CE for all possible Markov orders
for k = 0 : min(kMax, max(segLen) - 1)
    % A Markov process of order k is characterized by state sequences of
    % length k + 1; the current state plus the preceding k states.
    % collect all such subsequences of the data
    nSeq = sum(segLen - k);     % number of sequences of length k + 1
    seq = nan(nSeq, k + 1);     % sequences
    l = 0;
    for i = 1 : numel(segOff)   
        for j = 1 : segLen(i) - k
            l = l + 1;
            seq(l, :) = state(segOff(i) + (j : j + k))';
        end
    end
    % separate sequences into preceding k states and current state,
    preSeq = seq(:, 1 : k);
    curState = seq(:, k + 1);
    % and transform preceding sequence into integer
    [~, ~, preSeqID] = unique(preSeq, 'rows');
    % count transitions from preceding sequence to current state
    c = accumarray([preSeqID, curState], 1);
    
    % transition probabilities
    P = bsxfun(@rdivide, c, sum(c, 2));
    if k == 1
        P1 = P;
    end
    % minus 2 log likelihood
    ind = (c > 0);
    m2lL(k + 1) = - sum(c(ind) .* log(P(ind))) + 0;     % trick to avoid -0
    
    % conditional entropy in bits
    CE(k + 1) = m2lL(k + 1) / nSeq / log(2);

    if max(sum(c, 2)) == 1
        % If from each preceding sequence of states, there is only one
        % possible current state, the estimated Markov process is
        % deterministic. At this point, m2lL has decreased to 0 and will
        % remain there. We can therefore leave the rest of the values at
        % their initialization value 0.
        break
    end
end

% orders
k = (0 : kMax)';

% BIC penalty
penalty = nStates .^ k * (nStates - 1) * log(n);

% BIC
BIC = penalty + m2lL;

% BIC estimate of Markov order
[~, ind] = min(BIC);
kEst = k(ind);
% fprintf('estimated Markov order is %d\n', kEst)

% construct results table
results = table(k, m2lL, penalty, BIC, CE);

if nargout == 0
    disp(results)
    clear kEst results
end