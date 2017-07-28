function [seq, ess, LL] = exactInferenceWithLL_hmfa(seq, params, varargin)
%
% [seq, ess, LL] = exactInferenceWithLL_hmfa(seq, params,...)
%
% Extracts latent trajectories given HMFA model parameters.
%
% INPUTS:
%
% seq         - data structure, whose nth entry (corresponding to the nth
%               experimental trial) has fields
%                 y (yDim x T) -- neural data
%                 T (1 x 1)    -- number of timesteps
% params      - HMFA model parameters
%  
% OUTPUTS:
%
% seq         - training data structure with new fields
%                 state (1 x T)             -- factor analyzer state
%                                              at each time point
%                 x (xDim x T x nStates)    -- latent neural state
%                                              at each time point
%                 p (nStates x T)           -- defined RECURSIVELY at 
%                                              time t (1 <= t <= T) as 
%                                              the probability of the 
%                                              most probable sequence of 
%                                              length t-1 for each factor
%                                              analyzer hidden state
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
%               for each factor analyzer (default: 1e6)
% getLL       - logical that specifies whether to compute
%               data log likelihood (default: false)
% getSeq      - logical that specifies whether to compute
%               seq fields (default: false)
%
% @ 2015 Akinyinka Omigbodun    aomigbod@ucsd.edu

  condNumLim                    = 1e6;
  getLL                         = false;
  getSeq                        = false;
  
  assignopts(who, varargin);
  
  nStates                       = params.nStates;
  faType                        = params.faType;
  
  d                             = params.d;
  C                             = params.C;
  R                             = params.R;

  [yDim, xDim, ~]               = size(C);

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

  emission.mu                   =...
    repmat(d, [1 (faType(1) + (1 - faType(1))*nStates)]);
  emission.Sigma                = zeros([yDim yDim nStates]);
  for j=1:nStates
    jC                          = j*faType(2) + (1 - faType(2));
    jR                          = j*faType(3) + (1 - faType(3));
  % temp                        = diag(1./diag(R(:,:,jR)));
  % temp2                       = temp * C(:,:,jC);
  % emission.Sigma(:,:,j)       =...
  %   inv(temp - temp2 / (eye(xDim) + C(:,:,jC)' * temp2) * temp2');
    emission.Sigma(:,:,j)       = C(:,:,jC)*C(:,:,jC)' + R(:,:,jR);
    condNum                     = cond(emission.Sigma(:,:,j));
    if (condNum > condNumLim)
      error(['Covariance matrix of Gaussian distribution (state %d) ',...
             'has a large condition number (%d)'], j, condNum);
    end % if (condNum > condNumLim)
  end % for j=1:nStates

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
    if any(faType(2:3))
      if faType(3)
        invR                    = nan([yDim yDim nStates]);
      else
        invR                    = nan([yDim yDim]);
      end
      invRC                     = nan([yDim xDim nStates]);
      invM                      = nan([yDim yDim nStates]);
      betaFA                    = nan([xDim yDim nStates]);
    else % if ~any(faType(2:3))
      invRC                     = nan([yDim xDim]);
      invM                      = nan([yDim yDim]);
      betaFA                    = nan([xDim yDim]);
    end

    for n=1:N
      seq(n).x                  = zeros(xDim, seq(n).T, nStates);
    end % for n=1:N

    I                           = eye(xDim);
    for j=1:nStates
      jd                        = j*faType(1) + (1 - faType(1));
      jC                        = j*faType(2) + (1 - faType(2));
      jR                        = j*faType(3) + (1 - faType(3));
      jRC                       = max(jC,jR);

      invR(:,:,jR)              = diag(1./diag(R(:,:,jR)));
      invRC(:,:,jRC)            = invR(:,:,jR) * C(:,:,jC);
      invM(:,:,jRC)             =...
        invR(:,:,jR) - invRC(:,:,jRC) /...
                       (I + C(:,:,jC)' * invRC(:,:,jRC)) * invRC(:,:,jRC)';
      betaFA(:,:,jRC)           = C(:,:,jC)' * invM(:,:,jRC);

      for n=1:N
        seq(n).x(:,:,j)         =...
          betaFA(:,:,jRC)*bsxfun(@minus, seq(n).y, d(:,jd));      
      end % for n=1:N
    end % for j=1:nStates

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