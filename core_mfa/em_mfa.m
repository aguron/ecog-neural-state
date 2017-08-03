function [estParams, seq, LL, iterTime] = em_mfa(currentParams, seq, varargin)
%
% [estParams, seq, LL, iterTime] = em_mfa(currentParams, seq, ...)
%
% Fits MFA model parameters using expectation-maximization (EM) algorithm.
%
%   yDim: number of electrodes
%   xDim: state dimensionality
%
% INPUTS:
%
% currentParams - MFA model parameters at which EM algorithm is initialized
%                   faType (1 x 3)                  -- MFA factor
%                                                      analyzer(s)
%                   nMixComp (1 x 1)                -- number of mixture
%                                                      components
%                   Pi (1 x nMixComp)               -- mixture component
%                                                      priors
%                   d (yDim x nMixComp)             -- observation mean(s)
%                   C (yDim x xDim x nMixComp)      -- factor loadings
%                   R (yDim x yDim x nMixComp (1))	-- observation noise
%                                                      covariance(s)
% seq           - training data structure, whose nth entry (corresponding
%                 to the nth experimental trial) has fields
%                   trialId                     -- unique trial identifier
%                   segId                       -- segment identifier
%                                                  within trial
%                   T (1 x 1)                   -- number of timesteps in 
%                                                  segment
%                   y (yDim x T)                -- ECoG data
%
% OUTPUTS:
%
% estParams     - learned MFA model parameters returned by EM algorithm
%                   (same format as currentParams)
% seq           - training data structure with new fields
%                   mixComp (1 x T)             -- most probable factor
%                                                  analyzer mixture
%                                                  component at each
%                                                  time point
%                   x (xDim x T x nMixComp)     -- latent neural state at
%                                                  each time point
%                   p (nMixComp x T)            -- factor analyzer mixture
%                                                  component posterior
%                                                  probabilities at each
%                                                  time point
% LL            - data log likelihood after each EM iteration
% iterTime      - computation time for each EM iteration
%               
% OPTIONAL ARGUMENTS:
%
% emMaxIters    - number of EM iterations to run (default: 500)
% tolMFA        - stopping criterion for EM (default: 1e-2)
% minVarFrac    - fraction of overall data variance for each observed
%                 dimension to set as the private variance floor.  This 
%                 is used to combat Heywood cases, where ML parameter 
%                 learning returns one or more zero private variances.
%                 (default: 0.01)
%                 (See Martin & McDonald, Psychometrika, Dec 1975.)
% verbose       - logical that specifies whether to display status messages
%                 (default: false)
% freqLL        - data likelihood is computed every freqLL EM iterations. 
%                 freqLL = 1 means that data likelihood is computed every 
%                 iteration. (default: 1)
%
% Code adapted from mfa.m by Zoubin Ghahramani.
%
% @ 2016 Akinyinka Omigbodun    aomigbod@ucsd.edu

  emMaxIters                  	= 500;
  tolMFA                        = 1e-2;
  minVarFrac                    = 0.01;
  verbose                       = false;
  freqLL                        = 1;
  
  dbg                          	= false;

  assignopts(who, varargin);

  nMixComp                      = currentParams.nMixComp;
  faType                        = currentParams.faType;

  if ~isequal(faType,[1 1 1]) &&...
     ~isequal(faType,[1 1 0]) &&...
     ~isequal(faType,[1 0 0])
   fprintf('Does not support faType = [%d %d %d]\n',faType);
  end

  [yDim, xDim, ~]               = size(currentParams.C);

  if any(faType(2:3))
    if faType(3)
      invR                      = nan([yDim yDim nMixComp]);
    else % if ~faType(3)
      invR                      = nan([yDim yDim]);
    end
    invRC                       = nan([yDim xDim nMixComp]);
    invM                        = nan([yDim yDim nMixComp]);
    beta                        = nan([xDim yDim nMixComp]);
  else % if ~any(faType(2:3))
    invRC                       = nan([yDim xDim]);
    invM                        = nan([yDim yDim]);
    beta                        = nan([xDim yDim]);
  end
  YY_dd                         = nan([yDim yDim nMixComp]);
  YY_ddBeta                    	= nan([yDim xDim nMixComp]);
  EXX                          	= nan([xDim xDim nMixComp]);

  yAll                         	= [seq.y];
  lY                           	= size(yAll, 2);
  H                           	= nan(lY, nMixComp);	% E(w|y)
  logH                         	= nan(lY, nMixComp);	% log[E(w|y)]
  EX                           	= nan([xDim lY nMixComp]);
% const                        	= (2*pi)^(-yDim/2);
  logconst                    	= (-yDim/2)*log(2*pi);
  I                            	= eye(xDim);

  LL                           	= [];
  LLi                          	= -Inf;
  iterTime                    	= [];
  varFloor                    	= minVarFrac * diag(cov(yAll', 1));

  flag_converged               	= 0;
  
  % Loop once for each iteration of EM algorithm
  for i=1:max(emMaxIters,1)
    if verbose
      fprintf('\n');
    end % if verbose
    tic;

    if (emMaxIters)
      fprintf('EM iteration %3d of %d', i, emMaxIters);
    end % if (emMaxIters)
    if (rem(i, freqLL) == 0) || (i<=2)
      getLL                   	= true;
    else
      getLL                    	= false;
    end

    if ~isnan(LLi)
      LLold                    	= LLi;
    end % if ~isnan(LLi)

    % ==== E STEP =====
    d                         	= currentParams.d;
    C                          	= currentParams.C;
    R                          	= currentParams.R;
    Pi                         	= currentParams.Pi;
    for j=1:nMixComp
      jd                       	= j*faType(1) + (1 - faType(1));
      jC                       	= j*faType(2) + (1 - faType(2));
      jR                       	= j*faType(3) + (1 - faType(3));
      jRC                      	= max(jC,jR);

      invR(:,:,jR)             	= diag(1./diag(R(:,:,jR)));
      invRC(:,:,jRC)           	= invR(:,:,jR) * C(:,:,jC);
      invM(:,:,jRC)            	= invR(:,:,jR) - invRC(:,:,jRC) /...
                                  (I + C(:,:,jC)' * invRC(:,:,jRC)) *...
                                  invRC(:,:,jRC)';
      ldM                      	= sum(log(diag(chol(invM(:,:,jRC)))));
      S                        	= bsxfun(@minus,yAll,d(:,jd));
      S                       	= S';
      SinvM                    	= S * invM(:,:,jRC);
%     H(:,j)                   	= const*Pi(j)*...
%                                 exp(ldM - 0.5*sum(SinvM.*S, 2));
      logH(:,j)                	= logconst + log(Pi(j)) +...
                                  ldM - 0.5*sum(SinvM.*S, 2);
      EX(:,:,j)                	= (SinvM * C(:,:,jC))';
    end % for j=1:nMixComp

%   Hzero                      	= (sum(H,2)==0);
%   H(Hzero,:)                 	= eps(0);

%   H                          	= H + eps(0);    
%   LLi                        	= sum(log(sum(H,2)));
%   H                          	= normalize(H,2);
    
    [logH, scale]             	= normalizeLogspace(logH);
    LLi                       	= sum(scale);
    H                          	= exp(logH);

    Hsum                       	= sum(H);
    Hsumsum                    	= sum(Hsum);

    for j=1:nMixComp
      jd                       	= j*faType(1) + (1 - faType(1));
      jC                       	= j*faType(2) + (1 - faType(2));
      jR                       	= j*faType(3) + (1 - faType(3));
      jRC                      	= max(jC,jR);
      
      invR(:,:,jR)             	= diag(1./diag(R(:,:,jR)));
      invRC(:,:,jRC)           	= invR(:,:,jR) * C(:,:,jC);
      invM(:,:,jRC)            	= invR(:,:,jR) - invRC(:,:,jRC) /...
                                  (I + C(:,:,jC)' * invRC(:,:,jRC)) *...
                                  invRC(:,:,jRC)';
      S                        	= bsxfun(@dotprod,...
                                       	 bsxfun(@minus,yAll,d(:,jd)),...
                                       	 sqrt(H(:,j))');
      S                       	= S';
      YY_dd(:,:,j)             	= (S'*S)/Hsum(j);
      beta(:,:,jRC)            	= C(:,:,jC)' * invM(:,:,jRC);
      YY_ddBeta(:,:,j)        	= YY_dd(:,:,j) * beta(:,:,jRC)';
      EXX(:,:,j)              	= I - beta(:,:,jRC) * C(:,:,jC) +...
                                  beta(:,:,jRC) * YY_ddBeta(:,:,j);
    end % for j=1:nMixComp

    LL                        	= [LL LLi];

    if ~(emMaxIters)
      break
    end % if ~(emMaxIters)
    
    % Verify that likelihood is growing monotonically
    if (LLi < LLold)
      if (dbg)
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
           ((LLi-LLbase) < (1+tolMFA)*(LLold-LLbase))
      dispLL(false, getLL, LLi);
      flag_converged            = true;
      break
    else
      if (i<=2)
        LLbase                  = LLi;
      end % if (i<=2)
      previousParams            = currentParams;
    end

    % ==== M STEP =====
    if ~faType(2)
      if faType(1)
        currentParams.d        	= bsxfun(@rdivide, yAll*H, Hsum);
      end % if faType(1)
      currentParams.C         	=...
       sum(bsxfun(@times,YY_ddBeta,reshape(Hsum,1,1,[])), 3) /...
       sum(bsxfun(@times,EXX,reshape(Hsum,1,1,[])), 3);
    end % if ~faType(2)

    currentParams.R(:)         	= 0;
    for j=1:nMixComp
      T0                        = bsxfun(@dotprod, yAll, H(:,j)');
      T1                        = T0*[EX(:,:,j)' ones(lY, 1)];
      T2                        = [Hsum(j)*EXX(:,:,j) EX(:,:,j)*H(:,j);
                                   H(:,j)'*EX(:,:,j)' Hsum(j)];
      T3                        = T1 / T2;
      if faType(2)
        if faType(1)
          currentParams.d(:,j) 	= T3(:,xDim+1);
        end % if faType(1)
        currentParams.C(:,:,j) 	= T3(:,1:xDim);
      end % if faType(2)

      if faType(3)
        if currentParams.notes.RforceDiagonal
          diagR                	= diag(T0*yAll' - T3*T1')/Hsum(j);
                                % Set minimum private variance
          diagR               	= max(varFloor, diagR);
          currentParams.R(:,:,j)= diag(diagR);
        else
          currentParams.R(:,:,j)= (T0*yAll' - T3*T1')/Hsum(j);
                                % ensure symmetry
          currentParams.R(:,:,j)= symm(currentParams.R(:,:,j));
        end
      else % if ~faType(3)
        if currentParams.notes.RforceDiagonal
          diagR               	= diag(T0*yAll' - T3*T1')/Hsumsum;
                                % Set minimum private variance
          diagR                	= max(varFloor, diagR);
          currentParams.R      	= currentParams.R + diag(diagR);
        else
          currentParams.R      	= currentParams.R +...
                                  symm((T0*yAll' - T3*T1')/Hsumsum);
        end
      end
    end % for j=1:nMixComp

    currentParams.Pi          	= Hsum/Hsumsum;

    tEnd                       	= toc;
    iterTime                  	= [iterTime tEnd];
    
    % Display the most recent likelihood that was evaluated
    dispLL(verbose, getLL, LLi, tEnd);
  end % for i=1:max(emMaxIters,1)

  seq                          	= segmentByTrial(seq, H', 'p');
  seq                         	= segmentByTrial2(seq, EX, 'x');
  seq                         	=...
    segmentByTrial(seq, maxidx(H',[],1), 'mixComp');

  if (emMaxIters)
    fprintf('\n');
  end % if (emMaxIters)
  if flag_converged
    fprintf('Fitting has converged after %d EM iterations.\n', numel(LL));
  end % if flag_converged

  for j=1:size(currentParams.R,3)
    if any(diag(currentParams.R(:,:,j)) == varFloor)
      fprintf(['Warning: Private variance floor used for one or more ',...
               'observed dimensions (mixture component %d) in MFA.\n'], j);
    end % if any(diag(currentParams.R(:,:,j)) == varFloor)
  end % for j=1:size(currentParams.R,3)

  estParams                   	= currentParams;
end


% Display the most recent likelihood that was evaluated
function dispLL(verbose, getLL, LLi, tEnd)
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
end