function [estParams, seq, LL, iterTime] = em_mhmm(currentParams, seq, varargin)
%
% [estParams, seq, LL, iterTime] = em_mhmm(currentParams, seq, ...)
%
% Fits MHMM model parameters using expectation-maximization (EM) algorithm.
%
%   yDim: number of electrodes
%
% INPUTS:
%
% currentParams - HMM model parameters at which EM algorithm is initialized
%                   covType                         -- HMM covariance type
%                   nStates (1 x 1)                 -- number of HMM states
%                   pi (1 x nStates)                -- start probabilities
%                   piPrior (1 x nStates)           -- pseudo counts for HMM
%                                                      start probabilities
%                   trans (nStates x nStates)       -- transition matrix
%                   transPrior (nStates x nStates)  -- pseudo counts for HMM
%                                                      transition matrix
%                   d (yDim x nStates)              -- observation mean(s)
%                   R (yDim x yDim x nStates)       -- observation noise
%                                                      covariance(s)
% seq           - training data structure, whose nth entry (corresponding to
%                 the nth experimental trial) has fields
%                   trialId      -- unique trial identifier
%                   segId        -- segment identifier within trial
%                   T (1 x 1)    -- number of timesteps in segment
%                   y (yDim x T) -- ECoG data
%
% OUTPUTS:
%
% estParams     - learned HMM model parameters returned by EM algorithm
%                   (same format as currentParams)
% seq           - training data structure with new fields
%                 state (1 x T)             -- state at each time point
%                 mixComp (1 x 1)           -- most probable factor 
%                                              MHMM mixture component
%                 p (nStates x T)           -- state probabilities at each
%                                              time point (please see
%                                              EXACTINFERENCEWITHLL_MHMM
%                                              for details)
%                 P (1 x nMixComp)         	-- MHMM mixture component
%                                              posterior probabilities
% LL            - data log likelihood after each EM iteration
% iterTime      - computation time for each EM iteration
%               
% OPTIONAL ARGUMENTS:
%
% emMaxIters    - number of EM iterations to run (default: 500)
% tolHMM        - stopping criterion for EM (default: 1e-2)
% minVarFrac    - fraction of overall data variance for each observed dimension
%                 to set as the private variance floor.  This is used to combat
%                 Heywood cases, where ML parameter learning returns one or more
%                 zero private variances. (default: 0.01)
%                 (See Martin & McDonald, Psychometrika, Dec 1975.)
% verbose       - logical that specifies whether to display status messages
%                 (default: false)
% freqLL        - data likelihood is computed every freqLL EM iterations. 
%                 freqLL = 1 means that data likelihood is computed every 
%                 iteration. (default: 1)
%
% Code adapted from em.m by Byron Yu and John Cunningham.
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  emMaxIters                        = 500;
  tolHMM                            = 1e-2;
  verbose                           = false;
  freqLL                            = 1;
  outliers                          = [];

  extraOpts                        	= assignopts(who, varargin);

  nStates                           = currentParams(1).nStates;
  covType                          	= currentParams(1).covType;
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
      fprintf(['\nError: Data likelihood has decreased from %g to %g.\n',...
               'Please check component Gaussians'' positive-definiteness.\n'],...
              LLold, LLi);
      fprintf(['To continue running EM algorithm, execute USERCONTROL\n',...
               'in a separate MATLAB instance AND the SAME directory\n'...
               'where EM_MHMM is being executed before execution gets\n',...
               'here. When the execution stops, set dbg = true, and\n',...
               'execute DBCONT\n']);
      dbg                           = false;
      programcontrol
      if (dbg)
        continue
      else
        currentParams               = previousParams;
        flag_converged             	= true;
        break
      end % if (dbg)
    elseif exist('LLbase','var') &&...
           ((LLi-LLbase) < (1+tolHMM)*(LLold-LLbase))
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