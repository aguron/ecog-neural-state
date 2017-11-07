function [estParams, seq, LL, iterTime] = em_hmm(currentParams, seq, varargin)
%
% [estParams, seq, LL, iterTime] = em_hmm(currentParams, seq, ...)
%
% Fits HMM model parameters using expectation maximization
%
%   yDim: number of electrodes or channels
%
% INPUTS:
%
%
% currentParams - HMM model parameters with which EM algorithm is
%                 initialized in the fields
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
%                   d (yDim x nStates)                -- observation means
%                   R (yDim x yDim x nStates (or 1))	-- covariance(s)
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
% estParams     - learned HMM model parameters returned by EM algorithm
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
%                   state (1 x T)           -- HMM state at each time
%                                              point (from the Viterbi
%                                              path)
%                   p (nStates x T)       	-- hidden state probabilities
%                                              at each time point (please
%                                              see EXACTINFERENCEWITHLL_HMM
%                                              for details)
% LL            - data loglikelihood after each EM iteration
% iterTime      - computation time for each EM iteration
%
% OPTIONAL ARGUMENTS:
%
% emMaxIters    - number of EM iterations to run (default: 500)
% tolHMM        - stopping criterion for EM (default: 1e-2)
% verbose       - logical that specifies whether to display status messages
%                 (default: false)
% freqLL        - data loglikelihood is computed every freqLL EM
%                 iterations (default: 1)
%
% dbg           - set to true for EM debugging mode (default: false)
%
% Code adapted from em.m by Byron Yu and John Cunningham.
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  emMaxIters                        = 500;
  tolHMM                            = 1e-2;
  verbose                           = false;
  freqLL                            = 1;
  
  dbg                               = false;

  extraOpts                        	= assignopts(who, varargin);

  nStates                           = currentParams.nStates;
  covType                           = currentParams.covType;
  sharedCov                         = currentParams.sharedCov;
  
  yDim                              = size(currentParams.R, 1);
  
  YY_dd                             = nan([yDim yDim nStates]);

  LL                                = [];
  LLi                               = -Inf;
  iterTime                          = [];

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
    end % if ~isnan(LLi)

    % ==== E STEP =====
    [~, ess, LLi]                   =...
      exactInferenceWithLL_hmm(seq, currentParams,...
                               'getLL', getLL,...
                               extraOpts{:});
    LL                              = [LL LLi];

    % ==== M STEP =====
    currentParams.pi                =...
      normalize(ess.startCounts + currentParams.piPrior);
    currentParams.trans             =...
      normalize(ess.transCounts + currentParams.transPrior, 2);
    currentParams.d                 = ess.ybar;
    d                               = currentParams.d;
    for j=1:nStates
      S                             =...
        bsxfun(@dotprod,...
               bsxfun(@minus,[seq.y],d(:,j)), sqrt(ess.weights(:,j))');
      S                             = S';
      YY_dd(:,:,j)                  = (S'*S)/ess.wsum(j);
    end % for j=1:nStates

    if isequal(covType, 'full')
      currentParams.R             	= YY_dd;
    elseif isequal(covType, 'diagonal')
      currentParams.R             	= nddiag(nddiag(YY_dd));
    end
    if (sharedCov)
     currentParams.R                =...
      sum(bsxfun(@times,currentParams.R,reshape(ess.wsum,1,1,[])),3)/...
      sum(ess.wsum);
    end % if (sharedCov)

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
        terminate                 	= false;
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
    elseif (nStates == 1) ||...
           (exist('LLbase','var') &&...
            ((LLi-LLbase) < (1+tolHMM)*(LLold-LLbase)))
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
    exactInferenceWithLL_hmm(seq, currentParams,...
                             'getSeq', true,...
                             extraOpts{:});
  
  fprintf('\n');
  if flag_converged
   fprintf('Fitting has converged after %d EM iteration(s).\n', numel(LL));
  end % if flag_converged

  estParams                         = currentParams;
end