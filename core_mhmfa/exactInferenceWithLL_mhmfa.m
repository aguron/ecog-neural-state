function [seq, ess, LL] = exactInferenceWithLL_mhmfa(seq, params, varargin)
%
% [seq, ess, LL] = exactInferenceWithLL_mhmfa(seq, params,...)
%
% Extracts latent trajectories given MHMFA model parameters.
%
% INPUTS:
%
% seq         - data structure, whose nth entry (corresponding to the nth
%               experimental trial) has fields
%                 y (yDim x T) -- neural data
%                 T (1 x 1)    -- number of timesteps
% params      - MHMFA model parameters
%  
% OUTPUTS:
%
% seq         - training data structure with new fields
%                 state (1 x T)             -- factor analyzer state
%                                              at each time point
%                 mixComp (1 x 1)           -- most probable factor 
%                                              MHMFA mixture component
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
%                 P (1 x nMixComp)         	-- MHMFA mixture component
%                                              posterior probabilities
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
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  condNumLim                    = 1e6;
  getLL                         = false;
  getSeq                        = false;
  parallel                      = false;
  
  assignopts(who, varargin);
  
  nStates                       = params(1).nStates;
  faType                        = params(1).faType;
  
  nMixComp                     	= params(1).nMixComp;
  
  [yDim, xDim, ~]               = size(params(1).C);

  % Compute the expected sufficient statistics
  N                             = numel(seq);
  yAll                          = [seq.y];
  lY                            = size(yAll, 2);
  seqidx                        = cumsum([1, [seq.T]]); % keep track of
                                                        % where sequences
                                                        % start

  ess                           = struct('startCounts',cell(1,nMixComp),...
                                         'transCounts',cell(1,nMixComp),...
                                         'weights',cell(1,nMixComp));
  [ess(:).startCounts]          = deal(zeros(1, nStates));
  [ess(:).transCounts]          = deal(zeros(nStates));
  [ess(:).weights]              = deal(zeros(lY, nStates));
  ess(1).H                     	= nan(N, nMixComp);

  if (getLL)
    LL                          = 0;
  else % if (~getLL)
    LL                          = NaN;
  end

  emission                      = struct('mu',cell(1,nMixComp),...
                                         'Sigma',cell(1,nMixComp),...
                                         'd',cell(1,nMixComp),...
                                         'cpdType',cell(1,nMixComp),...
                                         'nstates',cell(1,nMixComp));
  for k=1:nMixComp
    emission(k).mu             	=...
      repmat(params(k).d, [1 (faType(1) + (1 - faType(1))*nStates)]);
    emission(k).Sigma         	= nan([yDim yDim nStates]);
    for j=1:nStates
      jC                       	= j*faType(2) + (1 - faType(2));
      jR                        = j*faType(3) + (1 - faType(3));
    % temp                      = diag(1./diag(R(:,:,jR)));
    % temp2                     = temp * C(:,:,jC);
    % emission.Sigma(:,:,j)     =...
    %   inv(temp - temp2 / (eye(xDim) + C(:,:,jC)' * temp2) * temp2');
      emission(k).Sigma(:,:,j)  =...
        params(k).C(:,:,jC)*params(k).C(:,:,jC)' + params(k).R(:,:,jR);
      condNum                   = cond(emission(k).Sigma(:,:,j));
      if (condNum > condNumLim)
        error(['Covariance matrix of Gaussian distribution',...
               '(state %d, mixComp %d) has a large',...
               'condition number (%d)'], j, k, condNum);
      end % if (condNum > condNumLim)
    end % for j=1:nStates
  end % for k=1:nMixComp

  [emission.d]                  = deal(yDim);
  [emission.cpdType]            = deal('condGauss');
  [emission.nstates]            = deal(nStates);

  logY                          = cell(1,nMixComp);
  for k=1:nMixComp
    logY{k}                    	= mkSoftEvidence(emission(k), yAll);
  end % for k=1:nMixComp
  [logY, scale]                 =...
    mufun({@normalizeLogspace, [true false]}, modifycells(logY,@(x)x'));
  Y                             = modifycells(logY,@(x)exp(x'));

  logp                          = nan(N,nMixComp);
  if (~parallel)
    for n=1:N
      ndx                     	= seqidx(n):seqidx(n+1)-1;

      for k=1:nMixComp
        Ynk                    	= Y{k}(:, ndx);
        [~, ~, ~, logp(n,k)]    =...
          hmmFwdBack2(params(k).pi, params(k).trans, Ynk);
      end % for k=1:nMixComp
    end % for n=1:N
  else % if (parallel)
    parfor k=1:nMixComp
    	for n=1:N
        ndx                    	= seqidx(n):seqidx(n+1)-1;
        Ynk                    	= Y{k}(:, ndx);
      
        [~, ~, ~, logp(n,k)]   	=...
          hmmFwdBack2(params(k).pi, params(k).trans, Ynk);
      end % for n=1:N
    end % parfor k=1:nMixComp
  end
  ess(1).H                     	= bsxfun(@plus,logp,log([params.Pi]));
  [ess(1).H, scale2]           	= normalizeLogspace(ess(1).H);
  ess(1).H                     	= exp(ess(1).H);
  p                             = ess(1).H;
  if (~parallel)
    for n=1:N
      ndx                     	= seqidx(n):seqidx(n+1)-1;

      for k=1:nMixComp
        Ynk                    	= Y{k}(:, ndx);
        [gamma, alpha, beta]    =...
          hmmFwdBack2(params(k).pi, params(k).trans, Ynk);
        ynSummed               	=...
          mhmmComputeTwoSliceSum(alpha, beta, params(k).trans, Ynk);
        ess(k).startCounts     	=...
          ess(k).startCounts + gamma(:, 1)' * p(n,k);
        ess(k).transCounts     	=...
          ess(k).transCounts + ynSummed * p(n,k);
        ess(k).weights(ndx, :) 	=...
          ess(k).weights(ndx, :) + gamma' * p(n,k);

      end % for k=1:nMixComp
    end % for n=1:N
  else % if (parallel)
    parfor k=1:nMixComp
    	for n=1:N
        ndx                    	= seqidx(n):seqidx(n+1)-1;
        Ynk                    	= Y{k}(:, ndx);
      
        [gamma, alpha, beta]   	=...
          hmmFwdBack2(params(k).pi, params(k).trans, Ynk);
        ynSummed               	=...
          mhmmComputeTwoSliceSum(alpha, beta, params(k).trans, Ynk);
        ess(k).startCounts     	=...
          ess(k).startCounts + gamma(:, 1)' * p(n,k);
        ess(k).transCounts     	=...
          ess(k).transCounts + ynSummed * p(n,k);
        ess(k).weights(ndx, :) 	=...
          ess(k).weights(ndx, :) + gamma' * p(n,k);
      end % for n=1:N
    end % parfor k=1:nMixComp
  end

  for k=1:nMixComp
    ess(k).Hsum                 = sum(ess(1).H(:,k));
    ess(k).wsum               	= sum(ess(k).weights, 1);
    ess(k).ybar                	=...yDim x nStates
      bsxfun(@rdivide, yAll*ess(k).weights, ess(k).wsum);
  end % for k=1:nMixComp

  if (getSeq)
    if any(faType(2:3))
      if faType(3)
        invR                    = nan([yDim yDim nStates nMixComp]);
      else % if ~faType(3)
        invR                    = nan([yDim yDim 1 nMixComp]);
      end
      invRC                     = nan([yDim xDim nStates nMixComp]);
      invM                      = nan([yDim yDim nStates nMixComp]);
      betaFA                    = nan([xDim yDim nStates nMixComp]);
    else
      invRC                     = nan([yDim xDim 1 nMixComp]);
      invM                      = nan([yDim yDim 1 nMixComp]);
      betaFA                    = nan([xDim yDim 1 nMixComp]);
    end

    for n=1:N
      seq(n).x                  = zeros(xDim, seq(n).T, nStates, nMixComp);
    end % for n=1:N

    I                           = eye(xDim);
    for j=1:nStates
      jd                        = j*faType(1) + (1 - faType(1));
      jC                        = j*faType(2) + (1 - faType(2));
      jR                        = j*faType(3) + (1 - faType(3));
      jRC                       = max(jC,jR);

      for k=1:nMixComp
        invR(:,:,jR,k)          = diag(1./diag(params(k).R(:,:,jR)));
        invRC(:,:,jRC,k)       	= invR(:,:,jR,k) * params(k).C(:,:,jC);
        invM(:,:,jRC,k)        	=...
        invR(:,:,jR,k) - invRC(:,:,jRC,k) /...
                         (I + params(k).C(:,:,jC)' * invRC(:,:,jRC,k)) *...
                         invRC(:,:,jRC,k)';
        betaFA(:,:,jRC,k)      	= params(k).C(:,:,jC)' * invM(:,:,jRC,k);

        for n=1:N
          seq(n).x(:,:,j,k)    	=...
            betaFA(:,:,jRC,k)*bsxfun(@minus, seq(n).y, params(k).d(:,jd));
        end % for n=1:N
      end % for k=1:nMixComp
    end % for j=1:nStates

    for n=1:N
      for k=1:nMixComp
        [seq(n).state(:,:,k), seq(n).p(:,:,k)]...
                                =hmmMap(struct('pi',params(k).pi,...
                                               'A', params(k).trans,...
                                               'emission', emission(k)),...
                                        seq(n).y);
      end % for k=1:nMixComp
      seq(n).mixComp            = maxidx(ess(1).H(n,:));
      seq(n).P                  = ess(1).H(n,:);
    end % for n=1:N
  end % if (getSeq)

  if (getLL)
    LL                          = LL + sum(scale) + sum(scale2);
    for k=1:nMixComp
      logprior                 	=...
        log(params(k).trans(:)+eps)'*(params(k).transPrior(:)) +...
        log(params(k).pi(:)+eps)'*(params(k).piPrior(:));
      LL                          = LL + logprior;
    end % for k=1:nMixComp
  end % if (getLL)
end