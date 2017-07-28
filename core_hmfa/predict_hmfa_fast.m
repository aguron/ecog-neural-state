function seq = predict_hmfa_fast(seq, params, varargin)
%
% seq = predict_hmfa_fast(seq, params, ...)
%
% Performs leave-electrode-out prediction for HMFA (assumes R is
% diagonal)
%
% INPUTS:
%
% seq         - test data structure
% params      - HMFA model parameters fit to training data
%
% OUTPUTS:
%
% seq      	- test data structure with new field
%               ycs (yDim x T)
%
% Code adapted from cosmoother_gpfa_viaOrth_fast.m by Byron Yu
% and John Cunningham.
%
% @ 2016 Akinyinka Omigbodun    aomigbod@ucsd.edu

  if ~params.notes.RforceDiagonal
    fprintf('ERROR: R must be diagonal to use predict_hmfa_fast.\n');
    return
  end % if ~params.notes.RforceDiagonal

  d                                 = params.d;
  C                                 = params.C;
  R                                 = params.R;
  
  pi                                = params.pi;
  trans                            	= params.trans;
  
  notes.RforceDiagonal              = params.notes.RforceDiagonal;

  [yDim, xDim, ~]                   = size(C);
  nStates                           = params.nStates;
  faType                            = params.faType;

  if any(faType(2:3))
    if faType(3)
      Rinv                          = nan([yDim yDim nStates]);
    else % if ~faType(3)
      Rinv                          = nan([yDim yDim]);
    end
    CRinv                           = nan([xDim yDim nStates]);
    CRinvC                          = nan([xDim xDim nStates]);
    B                               = nan([xDim xDim nStates]);
  else % if ~any(faType(2:3))
    Rinv                            = nan([yDim yDim]);
    CRinv                           = nan([xDim yDim]);
    CRinvC                          = nan([xDim xDim]);
    B                               = nan([xDim xDim]);
  end

  for j=1:nStates
    jC                              = j*faType(2) + (1 - faType(2));
    jR                              = j*faType(3) + (1 - faType(3));
    jRC                             = max(jC,jR);

    Rinv(:,:,jR)                    = diag(1./diag(R(:,:,jR)));
    CRinv(:,:,jRC)                  = C(:,:,jC)' * Rinv(:,:,jR);
    CRinvC(:,:,jRC)                 = CRinv(:,:,jRC) * C(:,:,jC);
    B(:,:,jRC)                      = inv(eye(xDim) + CRinvC(:,:,jRC));
    if ~any(faType(2:3))
      break
    end % if ~any(faType(2:3))
  end % for j=1:nStates

  Tall                              = [seq.T];

  for n=1:length(seq)
    T                               = Tall(n);
    ycs                             = nan(yDim, T);

    Yn                              = seq(n).y;
    parfor i=1:yDim
      % In computing the Viterbi path, eliminate the contribution
      % of electrode i
      mi                            = [1:(i-1) (i+1):yDim];

      [temp2, ~, ~]               	=...
        exactInferenceWithLL_hmfa(struct('y',Yn(mi,:), 'T',T),...
                                  struct('nStates',nStates,...
                                         'faType',faType,...
                                         'pi',pi, 'trans',trans,...
                                         'd',d(mi,:), 'C',C(mi,:,:),...
                                         'R',R(mi,mi,:)),...
                                  'getSeq', true,...
                                  varargin{:});
      state                         = temp2.state;

      if any(faType(2:3))
        temp                       	= mat2cell(B(:,:,state),...
                                               xDim, xDim, ones(1,T));
      else % if ~any(faType(2:3))
        temp                       	= cell(1,T);
        [temp{:}]                  	= deal(B);
      end

      % Taking advantage of the block diagonal matrix structure
      invM                          = blkdiag(temp{:});

      if faType(1)
                                    % yDim x T
        dif                         = Yn - d(:,state);
      else % if ~faType(1)
                                    % yDim x T
        dif                         = bsxfun(@minus, Yn, d);
      end

      if any(faType(2:3))
        CRinv_dif                   = nan(xDim, T);
        for t=1:T
          s                         = state(t);
          CRinv_dif(:,t)            = CRinv(:,:,s) * dif(:,t);
        end % for t=1:T
      else % if ~any(faType(2:3))
                                    % xDim x T
        CRinv_dif                   = CRinv * dif;
      end

      % Downdate invM to remove contribution of electrode i
      ci_invM                       = nan(T, xDim*T);
      ci_invM_ci                    = nan(T);
      idx                           = 1:xDim:(xDim*T + 1);
      if any(faType(2:3))
        ci                          = zeros([xDim nStates]);
        for j=1:nStates
          jC                        = j*faType(2) + (1 - faType(2));
          jR                        = j*faType(3) + (1 - faType(3));
          
          ci(:,j)                   = C(i,:,jC)' / sqrt(R(i,i,jR));
        end % for j=1:nStates
      else % if ~any(faType(2:3))
        ci                          = C(i,:)' / sqrt(R(i,i));
      end

      for t=1:T
        bIdx                        = idx(t):idx(t+1)-1;
        if any(faType(2:3))
          s                         = state(t);
          ci_invM(t,:)              = ci(:,s)' * invM(bIdx,:);
        else % if ~any(faType(2:3))
          ci_invM(t,:)              = ci' * invM(bIdx,:);
        end
      end % for t=1:T

      for t=1:T
        bIdx                        = idx(t):idx(t+1)-1;
        if any(faType(2:3))
          s                         = state(t);
          ci_invM_ci(:,t)           = ci_invM(:,bIdx) * ci(:,s);
        else % if ~any(faType(2:3))
          ci_invM_ci(:,t)           = ci_invM(:,bIdx) * ci;
        end
      end % for t=1:T
                                    % T x (xDim*T) 
      term                          = (ci_invM_ci - eye(T)) \ ci_invM;
                                    % (xDim*T) x (xDim*T)
      invM_mi                       = invM - ci_invM' * term;

      % Subtract out contribution of electrode i
      if any(faType(2:3))
        CRinvC_mi                   = nan([xDim xDim nStates]);
        for j=1:nStates
          CRinvC_mi(:,:,j)          = CRinvC(:,:,j) - ci(:,j) * ci(:,j)';
        end % for j=1:nStates
      else % if ~any(faType(2:3))
        CRinvC_mi                   = CRinvC - ci * ci';
      end

      if any(faType(2:3))
        term1Mat                    = nan([xDim T]);
        for t=1:T
          s                         = state(t);
          sC                        = s*faType(2) + (1 - faType(2));
          sR                        = s*faType(3) + (1 - faType(3));
          
          term1Mat(:,t)             =...
            CRinv_dif(:,t) - C(i,:,sC)' / R(i,i,sR) * dif(i,t);
        end % for t=1:T
                                    % (xDim*T) x 1
        term1Mat                    = term1Mat(:);
      else % if ~any(faType(2:3))
        term1Mat                    =... (xDim*T) x 1
          reshape(CRinv_dif - C(i,:)' / R(i,i) * dif(i,:), xDim*T, []);
      end

      % Compute blkProd = CRinvC_big * invM
      blkProd                       = nan(xDim*T, xDim*T);
      % idx                         = 1:xDim:(xDim*T + 1);
      for t=1:T
        bIdx                        = idx(t):idx(t+1)-1;
        if any(faType(2:3))
          s                         = state(t);
          blkProd(bIdx,:)           = CRinvC_mi(:,:,s) * invM_mi(bIdx,:);
        else % if ~any(faType(2:3))
          blkProd(bIdx,:)           = CRinvC_mi * invM_mi(bIdx,:);
        end
      end % for t=1:T
      blkProd                       = speye(xDim*T) - blkProd;
                                    % (xDim*T) x 1
      xsmMat                        = blkProd * term1Mat;
                                    % xDim x T
      xsmMat                        = reshape(xsmMat, xDim, T);

      if any(faType(1:2))
        temp3                       = nan(1,T);
        for t=1:T
          s                         = state(t);
          sd                        = s*faType(1) + (1 - faType(1));
          sC                        = s*faType(2) + (1 - faType(2));
          
          temp3(t)                  = C(i,:,sC) * xsmMat(:,t) + d(i,sd);
        end % for t=1:T
        ycs(i,:)                    = temp3;
      else % if ~any(faType(1:2))
        ycs(i,:)                    = C(i,:) * xsmMat + d(i);
      end
    end % parfor i=1:yDim
    seq(n).ycs                      = ycs;
    fprintf('Cross-validation complete for trial %d\r', n);
  end % for n=1:length(seq)
  fprintf('\n');
end