function [estParams, seq, LL, iterTime] = em_mhmm(currentParams, seq, varargin)
%
% [estParams, seq, LL, iterTime] = em_mhmm(currentParams, seq, ...)
%
% Fits MHMM model parameters using the Expectation Maximization (EM)
%	algorithm
%
%   yDim: number of electrodes or channels
%
% INPUTS:
%
% currentParams - MHMM model parameters with which EM algorithm is
%                 initialized in the fields (with the k-th entry of the
%                 struct corresponding to the k-th component HMM while
%                 faType, nMixComp, nStates, and notes are only specified
%                 in the 1st entry)
%                   nMixComp (1 x 1)                  -- number of
%                                                        component HMMs
%                   nStates (1 x 1)                   -- number of HMM
%                                                        states
%                   covType                           -- HMM component
%                                                        covariance type:
%                                                        'full' or
%                                                        'diagonal'
%                   sharedCov                         -- HMM components'
%                                                        covariance
%                                                        tied (true) or
%                                                        untied (false)
%                   pi (1 x nStates)                  -- start 
%                                                        probabilities
%                   piPrior (1 x nStates)             -- pseudo counts
%                                                        for start
%                                                        probabilities
%                   trans (nStates x nStates)         -- transition matrix
%                   transPrior (nStates x nStates)    -- pseudo counts for
%                                                        transition matrix
%                   Pi (1 x nMixComp)                 -- component HMFA
%                                                        priors
%                   d (yDim x nStates)                -- observation means
%                   R (yDim x yDim x nStates (or 1))	-- covariances
%
% seq           - training data structure, whose n-th entry (corresponding
%                 to the n-th experimental trial) has fields
%                   trialId         -- unique trial identifier
%                   segId           -- segment identifier within trial
%                   trialType       -- trial type index (Optional)
%                   fs              -- sampling frequency of ECoG data
%                   T (1 x 1)       -- number of timesteps in segment
%                   y (yDim x T)  	-- neural data
%
% OUTPUTS:
%
% estParams     - learned MHMM model parameters returned by EM algorithm
%                   (same format as currentParams)
% seq           - training data structure with fields
%                   trialId                 -- unique trial identifier
%                   segId                   -- segment identifier within 
%                                              trial
%                   trialType               -- trial type index (Optional)
%                   fs                      -- sampling frequency of ECoG
%                                              data
%                   T (1 x 1)               -- number of timesteps in
%                                              segment
%                   y (yDim x T)            -- neural data
%                   mixComp (1 x 1)        	-- most probable component HMM
%                   state (1 x T x         	-- HMM state at each time point
%                          nMixComp)           (from the Viterbi path)
%                   p (nStates x T          -- hidden state probabilities
%                      x nMixComp)             at each time point
%                                              (please see
%                                              EXACTINFERENCEWITHLL_MHMM
%                                              for details)
%                   P (1 x nMixComp)      	-- component HMM posterior
%                                              probabilities
% LL            - data loglikelihood after each EM iteration
% iterTime      - computation time for each EM iteration
%
% OPTIONAL ARGUMENTS:
%
% emMaxIters    - number of EM iterations to run (default: 500)
% tolMHMM     	- stopping criterion for EM (default: 1e-2)
% minVarFrac    - fraction of overall data variance for each observed
%                 dimension to set as the private variance floor.  This 
%                 is used to combat Heywood cases, where ML parameter 
%                 learning returns one or more zero private variances.
%                 (default: 0.01)
%                 (See Martin & McDonald, Psychometrika, Dec 1975.)
% verbose       - logical that specifies whether to display status messages
%                 (default: false)
% freqLL        - data loglikelihood is computed every freqLL EM 
%                 iterations (default: 1)
% outliers      - vector of indices of outlier trials in seq (default: [])
%
% dbg           - set to true for EM debugging mode (default: false)
%
% Code adapted from em.m by Byron Yu and John Cunningham.
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  emMaxIters                        = 500;
  tolMHMM                            = 1e-2;
  verbose                           = false;
  freqLL                            = 1;
  outliers                          = [];

  dbg                               = false;

  extraOpts                        	= assignopts(who, varargin);

  nStates                           = currentParams(1).nStates;
  covType                          	= currentParams(1).covType;
  sharedCov                         = currentParams(1).sharedCov;
  nMixComp                          = currentParams(1).nMixComp;
  
  yDim                              = size(currentParams(1).R, 1);

  YY_dd                             = nan([yDim yDim nStates nMixComp]);

  LL                                = [];
  LLi                               = -Inf;
  iterTime                          = [];
  inliers                           = setdiff(1:numel(seq),outliers);

  flag_converged                    = 0;

  % Loop once for each iteration of EM algorithm
  for i=1:emMaxIters
    if verbose
      fprintf('\n');
    end % if verbose
    tic;

    fprintf('EM iteration %3d of %d', i, emMaxIters);
    if (rem(i, freqLL) == 0) || (i<=2)
      getLL                         = true;
    else
      getLL                         = false;
    end

    if ~isnan(LLi)
      LLold                         = LLi;
    end

    % ==== E STEP =====
    [~, ess, LLi]                  	=...
      exactInferenceWithLL_mhmm(seq(inliers), currentParams,...
                                'getLL', getLL,...
                                extraOpts{:});
    LL                              = [LL LLi];

    % ==== M STEP =====
    for k=1:nMixComp
      currentParams(k).pi         	=...
        normalize(ess(k).startCounts + currentParams(k).piPrior);
      currentParams(k).trans      	=...
        normalize(ess(k).transCounts + currentParams(k).transPrior, 2);
    end % for k=1:nMixComp
    Hsum                            = num2cell(normalize([ess.Hsum]));
    [currentParams(:).Pi]           = deal(Hsum{:});

    [currentParams(:).d]           	= deal(ess.ybar);

    for j=1:nStates
      for k=1:nMixComp
        S                           =...
          bsxfun(@dotprod,...
                 bsxfun(@minus,[seq(inliers).y],ess(k).ybar(:,j)),...
                 sqrt(ess(k).weights(:,j))');
        S                          	= S';
        YY_dd(:,:,j,k)              = (S'*S)/ess(k).wsum(j);
      end % for k=1:nMixComp
    end % for j=1:nStates

    for k=1:nMixComp
      if isequal(covType, 'full')
        currentParams(k).R        	= YY_dd(:,:,:,k);
      elseif isequal(covType, 'diagonal')
        currentParams(k).R         	= nddiag(nddiag(YY_dd(:,:,:,k)));
      end

      if (sharedCov)
       currentParams(k).R          	=...
        sum(bsxfun(@times,currentParams(k).R,reshape(ess(k).wsum,1,1,[])),3)/...
        sum(ess(k).wsum);
      end % if (sharedCov)
    end % for k=1:nMixComp
    
    tEnd                            = toc;
    iterTime                        = [iterTime tEnd];

    % Display the most recent likelihood that was evaluated
    if verbose
      if getLL
        fprintf('       lik %g (%.1f sec)\n', LLi, tEnd);
      else
        fprintf('\n');
      end
    else
      if getLL
        fprintf('       lik %g\r', LLi);
      else
        fprintf('\r');
      end
    end

    % Verify that likelihood is growing monotonically
    if (LLi < LLold)
      if (dbg)
        % For debugging
        fprintf(['Error: Data loglikelihood has decreased from %g',...
                 ' to %g.\nPlease check component Gaussians'' ',...
                 'positive-definiteness.\n'],...
                 LLold, LLi);
        fprintf('To stop the EM algorithm, set terminate = true \n');
        terminate                  	= false;
        keyboard
        if (terminate)
          currentParams           	= previousParams;
          flag_converged           	= true;
          break
        end % if (terminate)
      else % if (~dbg)
        error(['Error: Data loglikelihood has decreased from %g',...
               ' to %g.\nPlease check component Gaussians'' ',...
               'positive-definiteness.\n'],...
              LLold, LLi);
      end
    elseif exist('LLbase','var') &&...
           ((LLi-LLbase) < (1+tolMHMM)*(LLold-LLbase))
      flag_converged                = true;
      break
    else
      if (i<=2)
        LLbase                     	= LLi;
      end % if (i<=2)
      previousParams               	= currentParams;
    end
  end % for i=1:emMaxIters

  [seq, ~, ~]                       =...
    exactInferenceWithLL_mhmm(seq, currentParams,...
                              'getSeq', true,...
                              extraOpts{:});

  fprintf('\n');
  if flag_converged
    fprintf('Fitting has converged after %d EM iterations.\n', numel(LL));
  end % if flag_converged

  estParams                         = currentParams;
end