function [seq, ess, LL] = exactInferenceWithLL_hmm(seq, params, varargin)
%
% [seq, ess, LL] = exactInferenceWithLL_hmm(seq, params,...)
%
% Infers latent variables given HMM model parameters
%
% INPUTS:
%
% seq         - data structure, whose n-th entry (corresponding to the n-th
%               experimental trial) has fields
%               	trialId         -- unique trial identifier
%                	segId           -- segment identifier within trial
%                	trialType       -- trial type index (Optional)
%                	fs              -- sampling frequency of ECoG data
%                	T (1 x 1)       -- number of timesteps in segment
%                	y (yDim x T)  	-- neural data
% params      - HMM model parameters (same format as currentParams
%               in EM_HMM)
%
% OUTPUTS:
%
% seq        	- data structure with fields
%                 trialId                 -- unique trial identifier
%                 segId                   -- segment identifier within 
%                                            trial
%                 trialType               -- trial type index (Optional)
%                 fs                      -- sampling frequency of ECoG
%                                            data
%                 T (1 x 1)               -- number of timesteps in
%                                            segment
%                 y (yDim x T)            -- neural data
%                 state (1 x T)           -- HMM state at each time
%                                            point (from the Viterbi path)
%                 p (nStates x T)       	-- defined RECURSIVELY at 
%                                            time t (1 <= t <= T) as 
%                                            the probability of the 
%                                            most probable sequence of 
%                                            length t-1 for each hidden
%                                            state at time t, given the
%                                            observations from 1 to t
% ess         - expected sufficient statistics structure with fields
%                 startCounts             -- for start probabilities
%                 transCounts             -- for transition matrix
%                 weights                 -- smoothed marginals
%                 wsum                    -- sum of smoothed marginals over
%                                            time
%                 ybar                    -- observation means estimate
% LL          - data loglikelihood
%
% OPTIONAL ARGUMENTS:
%
% getLL       - logical that specifies whether to compute
%               data loglikelihood (default: false)
% getSeq      - logical that specifies whether to compute
%               new seq fields in output (default: false)
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  getLL                         = false;
  getSeq                        = false;
  
  assignopts(who, varargin);
  
  nStates                       = params.nStates;

  d                             = params.d;
  R                             = params.R;
  
  yDim                          = size(R,1);

  % Compute the expected sufficient statistics
  N                             = numel(seq);
  yAll                          = [seq.y];
  lY                            = size(yAll, 2);
  seqidx                        = cumsum([1, [seq.T]]); % keep track of
                                                        % where sequences
                                                        % start

  ess.startCounts               = zeros(1, nStates);
  ess.transCounts               = zeros(nStates, nStates);
  ess.weights                   = zeros(lY, nStates);
  
  if (getLL)
    LL                          = 0;
  else % if (~getLL)
    LL                          = NaN;
  end
  A                             = params.trans;
  pi                            = params.pi;

  emission.mu                   = d;
  emission.Sigma                = R;
  if (params.sharedCov)
   emission.Sigma               = repmat(emission.Sigma,[1 1 nStates]);
  end % if (params.sharedCov)

  emission.d                    = yDim;
  emission.cpdType              = 'condGauss';
  emission.nstates              = nStates;
  logY                          = mkSoftEvidence(emission, yAll);
  [logY, scale]                 = normalizeLogspace(logY');
  Y                             = exp(logY');
  for n=1:N
    ndx                         = seqidx(n):seqidx(n+1)-1;
    Yn                          = Y(:, ndx);
    [gamma, alpha, beta, logp]  = hmmFwdBack2(pi, A, Yn);
    if (getLL)
      LL                        = LL + logp;
    end % if (getLL)
    ynSummed                    = hmmComputeTwoSliceSum(alpha, beta, A, Yn);

    ess.startCounts             = ess.startCounts + gamma(:, 1)';
    ess.transCounts             = ess.transCounts + ynSummed;
    ess.weights(ndx, :)         = ess.weights(ndx, :) + gamma';
  end % for n=1:N
  ess.wsum                      = sum(ess.weights, 1);
  ess.ybar                      =...yDim x nStates
    bsxfun(@rdivide, yAll*ess.weights, ess.wsum);

  if (getSeq)
    for n=1:N
      [seq(n).state, seq(n).p]  =...
        hmmMap(structure(pi,A,emission), seq(n).y);
    end % for n=1:N    
  end % if (getSeq)

  if (getLL)
    LL                          = LL + sum(scale);
    logprior                    =...
      log(A(:)+eps)'*(params.transPrior(:)) +...
      log(pi(:)+eps)'*(params.piPrior(:));
    LL                          = LL + logprior;
  end % if (getLL)
end