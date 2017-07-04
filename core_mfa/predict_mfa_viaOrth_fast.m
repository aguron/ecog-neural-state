function seq = predict_mfa_viaOrth_fast(seq, params, mList, varargin)
%
% seq = predict_mfa_viaOrth_fast(seq, params, mList, ...)
%
% Performs leave-electrode-out prediction for MFA.  This version takes 
% advantage of R being diagonal for computational savings.
%
% INPUTS:
%
% seq         - test data structure
% params      - MFA model parameters fit to training data
% mList       - number of top orthonormal latent coordinates to use for 
%               prediction (e.g., 1:5)
%
% OUTPUTS:
%
% seq         - test data structure with new fields ycsOrthXX, where XX
%               are elements of mList.  seq(n).ycsOrthXX has the same
%               dimensions as seq(n).y.
%
% Code adapted from cosmoother_gpfa_viaOrth_fast.m by Byron Yu
% and John Cunningham.
%
% @ 2016 Akinyinka Omigbodun    aomigbod@ucsd.edu

  if ~params.notes.RforceDiagonal
    fprintf(['ERROR: R must be diagonal to use',...
             'predict_mfa_viaOrth_fast.\n']);
    return
  end % if ~params.notes.RforceDiagonal

  d                                 = params.d;
  C                                 = params.C;
  R                                 = params.R;

  [yDim, xDim, ~]                   = size(C);
  nMixComp                          = params.nMixComp;
  faType                            = params.faType;

  if any(faType(2:3))
    if faType(2)
      TT                            = nan([xDim xDim nMixComp]);
      Corth                         = nan([yDim xDim nMixComp]);
    else
      TT                            = nan([xDim xDim]);
      Corth                         = nan([yDim xDim]);
    end % if faType(2)
    if faType(3)
      Rinv                          = nan([yDim yDim nMixComp]);
    else
      Rinv                          = nan([yDim yDim]);
    end % if faType(3)
    CRinv                           = nan([xDim yDim nMixComp]);
    CRinvC                          = nan([xDim xDim nMixComp]);
    B                               = nan([xDim xDim nMixComp]);
  else
    TT                              = nan([xDim xDim]);
    Corth                           = nan([yDim xDim]);
    Rinv                            = nan([yDim yDim]);
    CRinv                           = nan([xDim yDim]);
    CRinvC                          = nan([xDim xDim]);
    B                               = nan([xDim xDim]);
  end % if any(faType(2:3))

  for j=1:nMixComp
    jC                              = j*faType(2) + (1 - faType(2));
    jR                              = j*faType(3) + (1 - faType(3));
    jRC                             = max(jC,jR);

    Rinv(:,:,jR)                    = diag(1./diag(R(:,:,jR)));
    CRinv(:,:,jRC)                  = C(:,:,jC)' * Rinv(:,:,jR);
    CRinvC(:,:,jRC)                 = CRinv(:,:,jRC) * C(:,:,jC);
    B(:,:,jRC)                      = inv(eye(xDim) + CRinvC(:,:,jRC));
    [~, Corth(:,:,jC), TT(:,:,jC)]  =...
      orthogonalize(zeros(xDim, 1), C(:,:,jC));

    if ~any(faType(2:3))
      break
    end % if ~any(faType(2:3))
  end % for j=1:nMixComp

  Tall                              = [seq.T];

  for n=1:length(seq)
    for m=mList
      fn                            = sprintf('ycsOrth%02d', m);
      seq(n).(fn)                   = nan(yDim, seq(n).T);
    end % for m=mList

    T                               = Tall(n);
    
    if ~any(faType(2:3))
      temp                          = cell(1,T);
      [temp{:}]                     = deal(B);
    end % if ~any(faType(2:3))

    Yn                              = seq(n).y;
    for i=1:yDim
      % In computing the most probable component for each time point, 
      % eliminate the contribution of electrode i
      seq(n).y(i,:)                 = [];
      
      params.d                      = d;
      params.d(i,:)                 = [];
      
      params.C                      = C;
      params.C(i,:,:)               = [];
      
      params.R                      = R;
      params.R(i,:,:)               = [];
      params.R(:,i,:)               = [];

      [~, temp2]                    =...
        em_mfa(params, seq(n), 'emMaxIters', 0);
      mixComp                       = temp2.mixComp;

      if any(faType(2:3))
        temp                       	= mat2cell(B(:,:,mixComp),...
                                               xDim, xDim, ones(1,T));
      end % if any(faType(2:3))

      % Taking advantage of the block diagonal matrix structure
      invM                          = blkdiag(temp{:});

      % Restore contribution of electrode i
      seq(n).y                      = Yn;
      if faType(1)
                                    % yDim x T
        dif                         = seq(n).y - d(:,mixComp);
      else
                                    % yDim x T
        dif                         = bsxfun(@minus, seq(n).y, d);
      end % if faType(1)
      
      if any(faType(2:3))
        CRinv_dif                   = nan(xDim, T);
        for t=1:T
          s                         = mixComp(t);
          CRinv_dif(:,t)            = CRinv(:,:,s) * dif(:,t);
        end % for t=1:T
      else % if ~any(faType(2:3))
                                    % xDim x T
        CRinv_dif                   = CRinv * dif;
      end % if any(faType(2:3))

      % Downdate invM to remove contribution of electrode i
      ci_invM                       = nan(T, xDim*T);
      ci_invM_ci                    = nan(T);
      idx                           = 1:xDim:(xDim*T + 1);
      if any(faType(2:3))
        ci                          = nan([xDim nMixComp]);
        for j=1:nMixComp
          jC                        = j*faType(2) + (1 - faType(2));
          jR                        = j*faType(3) + (1 - faType(3));
          
          ci(:,j)                   = C(i,:,jC)' / sqrt(R(i,i,jR));
        end % for j=1:nMixComp
      else
        ci                          = C(i,:)' / sqrt(R(i,i));
      end % if any(faType(2:3))

      for t=1:T
        bIdx                        = idx(t):idx(t+1)-1;
        if any(faType(2:3))
          s                         = mixComp(t);
          ci_invM(t,:)              = ci(:,s)' * invM(bIdx,:);
        else
          ci_invM(t,:)              = ci' * invM(bIdx,:);
        end % if any(faType(2:3))
      end % for t=1:T

      for t=1:T
        bIdx                        = idx(t):idx(t+1)-1;
        if any(faType(2:3))
          s                         = mixComp(t);
          ci_invM_ci(:,t)           = ci_invM(:,bIdx) * ci(:,s);
        else
          ci_invM_ci(:,t)           = ci_invM(:,bIdx) * ci;
        end % if any(faType(2:3))
      end % for t=1:T
                                    % T x (xDim*T) 
      term                          = (ci_invM_ci - eye(T)) \ ci_invM;
                                    % (xDim*T) x (xDim*T)
      invM_mi                       = invM - ci_invM' * term;

      % Subtract out contribution of electrode i
      if any(faType(2:3))
        CRinvC_mi                   = nan([xDim xDim nMixComp]);
        for j=1:nMixComp
          CRinvC_mi(:,:,j)          = CRinvC(:,:,j) - ci(:,j) * ci(:,j)';
        end % for j=1:nMixComp
      else
        CRinvC_mi                   = CRinvC - ci * ci';
      end % if any(faType(2:3))

      if any(faType(2:3))
        term1Mat                    = nan([xDim T]);
        for t=1:T
          s                         = mixComp(t);
          sC                        = s*faType(2) + (1 - faType(2));
          sR                        = s*faType(3) + (1 - faType(3));
          
          term1Mat(:,t)             =...
            CRinv_dif(:,t) - C(i,:,sC)' / R(i,i,sR) * dif(i,t);
        end % for t=1:T
                                    % (xDim*T) x 1
        term1Mat                    = term1Mat(:);
      else
        term1Mat                    =...(xDim*T) x 1
          reshape(CRinv_dif - C(i,:)' / R(i,i) * dif(i,:), xDim*T, []);
      end % if any(faType(2:3))

      % Compute blkProd = CRinvC_big * invM
      blkProd                       = nan(xDim*T, xDim*T);
      % idx                         = 1:xDim:(xDim*T + 1);
      for t=1:T
        bIdx                        = idx(t):idx(t+1)-1;
        if any(faType(2:3))
          s                         = mixComp(t);
          blkProd(bIdx,:)           = CRinvC_mi(:,:,s) * invM_mi(bIdx,:);
        else
          blkProd(bIdx,:)           = CRinvC_mi * invM_mi(bIdx,:);
        end % if any(faType(2:3))
      end % for t=1:T
      blkProd                       = speye(xDim*T) - blkProd;
                                    % (xDim*T) x 1
      xsmMat                        = blkProd * term1Mat;
                                    % xDim x T
      xsmMat                        = reshape(xsmMat, xDim, T);
      if faType(2)
        xorth                       = nan([xDim T]);
        for t=1:T
          s                         = mixComp(t);
          xorth(:,t)                = TT(:,:,s) * xsmMat(:,t);
        end % for t=1:T
      else
                                    % xDim x T
        xorth                       = TT * xsmMat;
      end % if faType(2)

      for m=mList
        fn                          = sprintf('ycsOrth%02d', m);
        if any(faType(1:2))
          for t=1:T
            s                       = mixComp(t);
            sd                      = s*faType(1) + (1 - faType(1));
            sC                      = s*faType(2) + (1 - faType(2));

            seq(n).(fn)(i,t)        =...
              Corth(i,1:m,sC) * xorth(1:m,t) + d(i,sd);
          end % for t=1:T
        else
          seq(n).(fn)(i,:)          = Corth(i,1:m) * xorth(1:m,:) + d(i);
        end % if any(faType(1:2))
      end % for m=mList
    end % for i=1:yDim
    fprintf('Cross-validation complete for trial %d\r', n);
  end % for n=1:length(seq)
  fprintf('\n');
end