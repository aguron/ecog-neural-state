function [estParams, seq, LL, iterTime] = em_hmfa(currentParams, seq, varargin)
%
% [estParams, seq, LL, iterTime] = em_hmfa(currentParams, seq, ...)
%
% Fits HMFA model parameters using expectation-maximization (EM) algorithm.
%
%   yDim: number of electrodes
%   xDim: state dimensionality
%
% INPUTS:
%
% currentParams - HMFA model parameters at which EM algorithm is initialized
%                   faType (1 x 3)                  -- HMFA factor analyzer(s)
%                   nStates (1 x 1)                 -- number of HMFA states
%                   pi (1 x nStates)                -- start probabilities
%                   piPrior (1 x nStates)           -- pseudo counts for HMM
%                                                      start probabilities
%                   trans (nStates x nStates)       -- transition matrix
%                   transPrior (nStates x nStates)  -- pseudo counts for HMM
%                                                      transition matrix
%                   d (yDim x nStates)              -- observation mean(s)
%                   C (yDim x xDim x nStates)       -- factor loadings (s)
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
% estParams     - learned HMFA model parameters returned by EM algorithm
%                   (same format as currentParams)
% seq           - training data structure with new fields
%                   state (1 x T)           -- factor analyzer state 
%                                              at each time point
%                   x (xDim x T x nStates)  -- latent neural state
%                                              at each time point
%                   p (nStates x T)         -- factor analyzer state
%                                              probabilities at each
%                                              time point (please see
%                                              EXACTINFERENCEWITHLL_HMFA
%                                              for details)
% LL            - data log likelihood after each EM iteration
% iterTime      - computation time for each EM iteration
%               
% OPTIONAL ARGUMENTS:
%
% emMaxIters    - number of EM iterations to run (default: 500)
% tolHMFA       - stopping criterion for EM (default: 1e-2)
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
% @ 2015 Akinyinka Omigbodun    aomigbod@ucsd.edu

  emMaxIters                        = 500;
  tolHMFA                           = 1e-2;
  minVarFrac                        = 0.01;
  verbose                           = false;
  freqLL                            = 1;

  extraOpts                        	= assignopts(who, varargin);

  nStates                           = currentParams.nStates;
  faType                            = currentParams.faType;

  if ~isequal(faType,[1 1 1]) &&...
     ~isequal(faType,[1 1 0]) &&...
     ~isequal(faType,[1 0 0])
   fprintf('Does not support faType = [%d %d %d]\n',faType);
  end

  [yDim, xDim, ~]                   = size(currentParams.C);

  if any(faType(2:3))
    if faType(3)
      invR                          = nan([yDim yDim nStates]);
    else
      invR                          = nan([yDim yDim]);
    end
    invRC                           = nan([yDim xDim nStates]);
    invM                            = nan([yDim yDim nStates]);
    beta                            = nan([xDim yDim nStates]);
  else % if ~any(faType(2:3))
    invRC                           = nan([yDim xDim]);
    invM                            = nan([yDim yDim]);
    beta                            = nan([xDim yDim]);
  end
  YY_dd                             = nan([yDim yDim nStates]);
  YY_ddBeta                         = nan([yDim xDim nStates]);
  EXX                               = nan([xDim xDim nStates]);

  LL                                = [];
  LLi                               = -Inf;
  iterTime                          = [];
  varFloor                          = minVarFrac * diag(cov([seq.y]', 1));

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
    [~, ess, ~]                     =...
      exactInferenceWithLL_hmfa(seq, currentParams, extraOpts{:});

    % ==== CM STEP 1 =====
    currentParams.pi                =...
      normalize(ess.startCounts + currentParams.piPrior);
    currentParams.trans             =...
      normalize(ess.transCounts + currentParams.transPrior, 2);
    if faType(1)
      currentParams.d               = ess.ybar;
    end % if faType(1)
    d                               = currentParams.d;

    % ==== E STEP =====
    [~, ess, LLi]                   =...
      exactInferenceWithLL_hmfa(seq, currentParams,...
                                'getLL', getLL,...
                                extraOpts{:});
    LL                              = [LL LLi];

    % ==== CM STEP 2 =====
    for j=1:nStates
      jd                            = j*faType(1) + (1 - faType(1));
      S                             =...
        bsxfun(@dotprod,...
               bsxfun(@minus,[seq.y],d(:,jd)), sqrt(ess.weights(:,j))');
      S                             = S';
      YY_dd(:,:,j)                  = (S'*S)/ess.wsum(j);
    end % for j=1:nStates
    I                               = eye(xDim);
    C                               = currentParams.C; % temporary
    R                               = currentParams.R; % temporary
    for j=1:nStates
      jC                            = j*faType(2) + (1 - faType(2));
      jR                            = j*faType(3) + (1 - faType(3));
      jRC                           = max(jC,jR);

      invR(:,:,jR)                  = diag(1./diag(R(:,:,jR)));
      invRC(:,:,jRC)                = invR(:,:,jR) * C(:,:,jC);
      invM(:,:,jRC)                 =...
        invR(:,:,jR) - invRC(:,:,jRC) /...
                     	 (I + C(:,:,jC)' * invRC(:,:,jRC)) *...
                       invRC(:,:,jRC)';
      beta(:,:,jRC)                 = C(:,:,jC)' * invM(:,:,jRC);
      YY_ddBeta(:,:,j)              = YY_dd(:,:,j) * beta(:,:,jRC)';
      EXX(:,:,j)                    = I - beta(:,:,jRC) * C(:,:,jC) +...
                                      beta(:,:,jRC) * YY_ddBeta(:,:,j);
    end % for j=1:nStates

    if faType(2)
      for j=1:nStates
        currentParams.C(:,:,j)    	= YY_ddBeta(:,:,j) / EXX(:,:,j);
      end % for j=1:nStates
    else % if ~faType(2)
      currentParams.C               =...
       sum(bsxfun(@times,YY_ddBeta,reshape(ess.wsum,1,1,[])), 3) /...
       sum(bsxfun(@times,EXX,reshape(ess.wsum,1,1,[])), 3);
    end

    C                               = currentParams.C; % temporary
    if faType(3)
      for j=1:nStates
        jC                         	= j*faType(2) + (1 - faType(2));
        if currentParams.notes.RforceDiagonal
          diagR                     =...
            diag(YY_dd(:,:,j)) - sum(YY_ddBeta(:,:,j) .* C(:,:,jC), 2);
          % Set minimum private variance
          diagR                     = max(varFloor, diagR);
          currentParams.R(:,:,j)    = diag(diagR);
        else
          currentParams.R(:,:,j)    =...
            YY_dd(:,:,j) - C(:,:,jC) * YY_ddBeta(:,:,j)';
          % ensure symmetry
          currentParams.R(:,:,j)    = symm(currentParams.R(:,:,j));
        end
      end % for j=1:nStates
    else % if ~faType(3)
      currentParams.R(:)            = 0;
      for j=1:nStates
        if currentParams.notes.RforceDiagonal
          diagR                     =...
            ess.wsum(j)*(diag(YY_dd(:,:,j)) -...
                         sum(YY_ddBeta(:,:,j) .* C(:,:,jC), 2));
          % Set minimum private variance
          diagR                     = max(varFloor, diagR);
          currentParams.R           = currentParams.R + diag(diagR);
        else
          % ensure symmetry
          currentParams.R           = currentParams.R +...
            ess.wsum(j)*symm(YY_dd(:,:,j) - C(:,:,jC) * YY_ddBeta(:,:,j)');
        end
      end % for j=1:nStates
      currentParams.R               = currentParams.R/sum(ess.wsum);
    end % if faType(3)

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
      fprintf('\nError: Data likelihood has decreased from %g to %g.\n',...
              LLold, LLi);
      fprintf('To continue EM algorithm, set dbg = true\n');
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
           ((LLi-LLbase) < (1+tolHMFA)*(LLold-LLbase))
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
    exactInferenceWithLL_hmfa(seq, currentParams,...
                              'getSeq', true,...
                              extraOpts{:});
  
  fprintf('\n');
  if flag_converged
    fprintf('Fitting has converged after %d EM iterations.\n', numel(LL));
  end % if flag_converged

  for j=1:size(currentParams.R,3)
    if any(diag(currentParams.R(:,:,j)) == varFloor)
      fprintf(['Warning: Private variance floor used for one ',...
               'or more observed dimensions (state = %d) in HMFA.\n'], j);
    end % if any(diag(currentParams.R(:,:,j)) == varFloor)
  end % for j=1:size(currentParams.R,3)

  estParams                         = currentParams;
end