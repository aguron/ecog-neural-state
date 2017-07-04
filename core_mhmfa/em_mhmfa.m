function [estParams, seq, LL, iterTime] = em_mhmfa(currentParams, seq, varargin)
%
% [estParams, seq, LL, iterTime] = em_mhmfa(currentParams, seq, ...)
%
% Fits MHMFA model parameters using expectation-maximization (EM) algorithm.
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
%                 state (1 x T)             -- factor analyzer state
%                                              at each time point
%                 mixComp (1 x 1)           -- most probable factor 
%                                              MHMFA mixture component
%                 x (xDim x T x nStates)    -- latent neural state
%                                              at each time point
%                 p (nStates x T)           -- factor analyzer state
%                                              probabilities at each
%                                              time point (please see
%                                              EXACTINFERENCEWITHLL_MHMFA
%                                              for details)
%                 P (1 x nMixComp)         	-- MHMFA mixture component
%                                              posterior probabilities
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
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  emMaxIters                        = 500;
  tolHMFA                           = 1e-2;
  minVarFrac                        = 0.01;
  verbose                           = false;
  freqLL                            = 1;
  outliers                          = [];

  extraOpts                        	= assignopts(who, varargin);

  nStates                           = currentParams(1).nStates;
  faType                            = currentParams(1).faType;
  nMixComp                          = currentParams(1).nMixComp;
  
  if ~isequal(faType,[1 1 1]) &&...
     ~isequal(faType,[1 1 0])
    fprintf('Does not support faType = [%d %d %d]\n',faType);
  end

  [yDim, xDim, ~]                   = size(currentParams(1).C);

  if any(faType(2:3))
    if faType(3)
      invR                          = nan([yDim yDim nStates nMixComp]);
    else
      invR                          = nan([yDim yDim 1 nMixComp]);
    end
    invRC                           = nan([yDim xDim nStates nMixComp]);
    invM                            = nan([yDim yDim nStates nMixComp]);
    beta                            = nan([xDim yDim nStates nMixComp]);
  else
    invRC                           = nan([yDim xDim 1 nMixComp]);
    invM                            = nan([yDim yDim 1 nMixComp]);
    beta                            = nan([xDim yDim 1 nMixComp]);
  end
  YY_dd                             = nan([yDim yDim nStates nMixComp]);
  YY_ddBeta                         = nan([yDim xDim nStates nMixComp]);
  EXX                               = nan([xDim xDim nStates nMixComp]);

  LL                                = [];
  LLi                               = -Inf;
  iterTime                          = [];
  inliers                           = setdiff(1:numel(seq),outliers);
  varFloor                          =...
    minVarFrac * diag(cov([seq(inliers).y]', 1));

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
    [~, ess, ~]                     =...
      exactInferenceWithLL_mhmfa(seq(inliers),currentParams,extraOpts{:});

    % ==== CM STEP 1 =====
    for k=1:nMixComp
      currentParams(k).pi         	=...
        normalize(ess(k).startCounts + currentParams(k).piPrior);
      currentParams(k).trans      	=...
        normalize(ess(k).transCounts + currentParams(k).transPrior, 2);
    end % for k=1:nMixComp
    Hsum                            = num2cell(normalize([ess.Hsum]));
    [currentParams(:).Pi]           = deal(Hsum{:});

    [currentParams(:).d]           	= deal(ess.ybar);

    % ==== E STEP =====
    [~, ess, LLi]                   =...
      exactInferenceWithLL_mhmfa(seq(inliers), currentParams,...
                                 'getLL', getLL,...
                                 extraOpts{:});
    LL                              = [LL LLi];

    % ==== CM STEP 2 =====
    for j=1:nStates
      jd                            = j*faType(1) + (1 - faType(1));
      for k=1:nMixComp
        S                           =...
          bsxfun(@dotprod,...
                 bsxfun(@minus,[seq(inliers).y],ess(k).ybar(:,jd)),...
                 sqrt(ess(k).weights(:,j))');
        S                          	= S';
        YY_dd(:,:,j,k)             	= (S'*S)/ess(k).wsum(j);
      end % for k=1:nMixComp
    end % for j=1:nStates
    I                               = eye(xDim);

    for j=1:nStates
      jC                            = j*faType(2) + (1 - faType(2));
      jR                            = j*faType(3) + (1 - faType(3));
      jRC                           = max(jC,jR);
      for k=1:nMixComp
        invR(:,:,jR,k)             	=...
          diag(1./diag(currentParams(k).R(:,:,jR)));
        invRC(:,:,jRC,k)           	=...
          invR(:,:,jR,k) * currentParams(k).C(:,:,jC);
        invM(:,:,jRC,k)           	=...
          invR(:,:,jR,k) -...
            invRC(:,:,jRC,k) /...
            (I + currentParams(k).C(:,:,jC)' * invRC(:,:,jRC,k)) *...
            invRC(:,:,jRC,k)';
        beta(:,:,jRC,k)           	=...
          currentParams(k).C(:,:,jC)' * invM(:,:,jRC,k);
          YY_ddBeta(:,:,j,k)       	= YY_dd(:,:,j,k) * beta(:,:,jRC,k)';
        EXX(:,:,j,k)               	=...
          I - beta(:,:,jRC,k) * currentParams(k).C(:,:,jC) +...
          beta(:,:,jRC,k) * YY_ddBeta(:,:,j,k);
      end % for k=1:nMixComp
    end % for j=1:nStates

    for j=1:nStates
      for k=1:nMixComp
        currentParams(k).C(:,:,j)   = YY_ddBeta(:,:,j,k)/EXX(:,:,j,k);
      end % for k=1:nMixComp
    end % for j=1:nStates

    if faType(3)
      for j=1:nStates
        for k=1:nMixComp
          if currentParams(1).notes.RforceDiagonal
           diagR                  	=...
            diag(YY_dd(:,:,j,k)) -...
            sum(YY_ddBeta(:,:,j,k) .* currentParams(k).C(:,:,j), 2);

           % Set minimum private variance
           diagR                    = max(varFloor, diagR);
           currentParams(k).R(:,:,j)= diag(diagR);
          else % if ~currentParams(1).notes.RforceDiagonal
           currentParams(k).R(:,:,j)=...
          	YY_dd(:,:,j,k) -...
            currentParams(k).C(:,:,j) * YY_ddBeta(:,:,j,k)';
          
           % ensure symmetry
           currentParams(k).R(:,:,j)= symm(currentParams(k).R(:,:,j));
          end
        end % for k=1:nMixComp
      end %  for j=1:nStates
    else % if ~faType(3)
      for k=1:nMixComp
        currentParams(k).R(:)       = 0;    
        for j=1:nStates
          if currentParams(1).notes.RforceDiagonal
            diagR                  	=...
             ess(k).wsum(j)*...
             (diag(YY_dd(:,:,j,k)) -...
              sum(YY_ddBeta(:,:,j,k) .* currentParams(k).C(:,:,j), 2));

            % Set minimum private variance
            diagR                  	= max(varFloor, diagR);
            currentParams(k).R     	= currentParams(k).R + diag(diagR);
          else % if ~currentParams(1).notes.RforceDiagonal
            currentParams(k).R     	= currentParams(k).R +...
              ess(k).wsum(j)*...
              symm(YY_dd(:,:,j,k) -...
                   currentParams(k).C(:,:,j) * YY_ddBeta(:,:,j,k)');
          end
        end % for j=1:nStates
        currentParams(k).R         	= currentParams(k).R/sum(ess(k).wsum);
      end % for k=1:nMixComp
    end

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
    exactInferenceWithLL_mhmfa(seq, currentParams,...
                               'getSeq', true,...
                               extraOpts{:});
  
  fprintf('\n');
  if flag_converged
    fprintf('Fitting has converged after %d EM iterations.\n', numel(LL));
  end % if flag_converged

  for j=1:size(currentParams(1).R,3)
   for k=1:nMixComp
    if any(diag(currentParams(k).R(:,:,j)) == varFloor)
      fprintf(['Warning: Private variance floor used for',...
               ' one or more observed dimensions (state %d,',...
               ' mixComp %d) in MHMFA.\n'], j, k);
    end % if any(diag(currentParams(k).R(:,:,j)) == varFloor)
   end % for k=1:nMixComp
  end % for j=1:size(currentParams(1).R,3)

  estParams                         = currentParams;
end