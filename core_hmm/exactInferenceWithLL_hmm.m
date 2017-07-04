function [seq, ess, LL] = exactInferenceWithLL_hmm(seq, params, varargin)
%
% [seq, ess, LL] = exactInferenceWithLL_hmm(seq, params,...)
%
% Extracts latent trajectories given HMM model parameters.
%
% INPUTS:
%
% seq         - data structure, whose nth entry (corresponding to the nth
%               experimental trial) has fields
%                 y (yDim x T) -- neural data
%                 T (1 x 1)    -- number of timesteps
% params      - HMM model parameters
%  
% OUTPUTS:
%
% seq         - training data structure with new fields
%                 state (1 x T)             -- hidden Markov state
%                                              at each time point
%                 p (nStates x T)           -- defined RECURSIVELY at 
%                                              time t (1 <= t <= T) as 
%                                              the probability of the 
%                                              most probable sequence
%                                              of length t-1 for each
%                                              hidden Markov state
%                                              at time t, given the
%                                              observations from 1 to t
% ess         - expected sufficient statistics structure with fields
%                 startCounts
%                 transCounts
%                 weights
%                 wsum
%                 ybar
% LL          - data log likelihood
%
% OPTIONAL ARGUMENTS:
%
% condNumLim  - upper limit of condition number of covariance
%               for each hidden Markov state (default: 1e6)
% getLL       - logical that specifies whether to compute
%               data log likelihood (default: false)
% getSeq      - logical that specifies whether to compute
%               seq fields (default: false)
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  condNumLim                    = 1e6;
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

  for j=1:max(1,nStates*~params.sharedCov)
    condNum                     = cond(emission.Sigma(:,:,j));
    if (condNum > condNumLim)
      error(['Covariance matrix of Gaussian distribution (state %d) ',...
             'has a large condition number (%d)'], j, condNum);
    end % if (condNum > condNumLim)
  end % for j=1:max(1,nStates*~params.sharedCov)
  
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