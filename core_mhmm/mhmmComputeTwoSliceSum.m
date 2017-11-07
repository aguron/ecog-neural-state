function xiSummed = mhmmComputeTwoSliceSum(alpha, beta, A, B)
%% Compute the sum of the two-slice distributions over hidden states for
%   a component HMM of an MHMM model, accounting for the possibility that
%   the probability of a sequence being generated (in its entirety) by
%   the component HMM is 0
%
% Let K be the number of hidden states, and T be the number of time steps.
% Let S(t) denote the hidden state at time t, and y(t) be the (not
% necessarily scalar) observation at time t. 
%
%% INPUTS:
% 
% alpha, and beta are computed using HMMFWDBACK2, A is the state
% transition matrix, whose *rows* sum to one, and B is the soft evidence. 
% 
% alpha(j, t)      = p( S(t) = j  | y(1:t)    )   (KxT) 
% beta (j, t) propto p( y(t+1:T)  | S(t)   = j)   (KxT)
% A    (i, j)      = p( S(t) = j  | S(t-1) = i)   (KxK) 
% B    (j, t)      = p( y(t)      | S(t)   = j)   (KxT)
% 
%% OUTPUT: 
%
% xiSummed(i, j) = sum_t=2:T p(S(t) = i, S(t+1) = j | y(1:T)), t=2:T	(KxK)
% The output constitutes the expected sufficient statistics for the 
% transition matrix, for a given observation sequence. 
%%
%
% Code adapted from hmmComputeTwoSliceSum.m in pmtk3
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

[K, T]        = size(B);
xiSummed      = zeros(K, K);
for t = T-1:-1:1
    b         = beta(:,t+1) .* B(:,t+1);
    xit      	= A .* (alpha(:,t) * b');
    if (sum(xit(:)) == 0)
    	continue
    end % if (sum(xit(:)) == 0)
    xiSummed  = xiSummed + xit./sum(xit(:));
end
end
