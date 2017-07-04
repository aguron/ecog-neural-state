function seq = predict_hmm(seq, params, varargin)
%
% seq = predict_hmm(seq, params, ...)
%
% Performs leave-electrode-out prediction for HMM.
%
% INPUTS:
%
% seq       	- test data structure
% params    	- HMM model fit to training data
%
% OUTPUTS:
%
% seq         - test data structure with new fields
%                 ycs (yDim x T)
%                 vcs (yDim x T)
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  d                                 = params.d;
  R                                 = params.R;

  pi                                = params.pi;
  trans                            	= params.trans;

  yDim                              = size(R,1);
  nStates                           = params.nStates;
  covType                           = params.covType;
  sharedCov                         = params.sharedCov;

  Tall                              = [seq.T];

  for n=1:numel(seq)
    T                               = Tall(n);
    ycs                             = nan(yDim, T);
    vcs                             = nan(yDim, T);
    Yn                              = seq(n).y;
    parfor i=1:yDim
      % In computing the Viterbi path, eliminate the contribution
      % of electrode i
      mi                            = [1:(i-1) (i+1):yDim];

      [temp2, ~, ~]                 =...
        exactInferenceWithLL_hmm(struct('y',Yn(mi,:), 'T',T),...
                                 struct('nStates',nStates,...
                                        'pi',pi, 'trans',trans,...
                                        'd',d(mi,:), 'R',R(mi,mi,:),...
                                        'sharedCov', sharedCov),...
                                 'getSeq', true,...
                                 varargin{:});
      state                         = temp2.state;
      
      % Taking advantage of block diagonal matrix structure
      invRmimi                    	= nan(yDim-1,yDim-1,nStates);
      for j=1:max(1,nStates*~sharedCov)
        if strcmp(covType, 'diagonal')
         invRmimi(:,:,j)          	= diag(diag(1./R(mi,mi,min(j,end))));
        elseif strcmp(covType, 'full')
         invRmimi(:,:,j)          	= inv(R(mi,mi,min(j,end)));
        end
      end % for j=1:max(1,nStates*~sharedCov)
      temp                          = mat2cell(invRmimi(:,:,min(state,end)),...
                                               yDim-1, yDim-1, ones(1,T));
                                    % (yDim-1)*T x (yDim-1)*T
      invBmimi                      = blkdiag(temp{:});

      if strcmp(covType, 'full')
        temp                        = mat2cell(R(i,mi,min(state,end)),...
                                               1, yDim-1, ones(1,T));
                                    % T x (yDim-1)*T
        Bimi                        = blkdiag(temp{:});
      end % if strcmp(covType, 'full')

      if strcmp(covType, 'full')
        temp                        = mat2cell(R(mi,i,min(state,end)),...
                                               yDim-1, 1, ones(1,T));
                                    % (yDim-1)*T x T
        Bmii                        = blkdiag(temp{:});
      end % if strcmp(covType, 'full')

                                    % T x T
      Bii                           = diag(matflat(R(i,i,min(state,end))));

      if strcmp(covType, 'diagonal')
        % 1 x T
        ycs(i,:)                    = d(i,state);
      elseif strcmp(covType, 'full')
                                    % yDim-1 x T
        difmi                       = Yn(mi,:) - d(mi,state);

        % 1 x T
        ycs(i,:)                    =...
          d(i,state) + (Bimi*invBmimi*difmi(:))';
      end
      
      if strcmp(covType, 'diagonal')
        % 1 x T
        vcs(i,:)                    = diag(Bii)';
      elseif strcmp(covType, 'full')
        % 1 x T
        vcs(i,:)                    = diag(Bii - Bimi*invBmimi*Bmii)';
      end
    end % parfor i=1:yDim
    seq(n).ycs                      = ycs;
    seq(n).vcs                      = vcs;
    fprintf('Cross-validation complete for trial %d\r', n);
  end % for n=1:numel(seq)
  fprintf('\n');
end