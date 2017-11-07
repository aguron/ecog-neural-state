function seq = predict_hmm(seq, params, varargin)
%
% seq = predict_hmm(seq, params, ...)
%
% Performs leave-one-channel-out prediction for HMM
%
% INPUTS:
%
% seq         - data structure (same format as seqTrain in HMMENGINE)
% params      - HMM model parameters (same format as currentParams
%               in EM_HMM)
%
% OUTPUTS:
%
% seq        	- data structure with fields
%                 trialId                 -- unique trial identifier
%                 segId                   -- segment identifier within 
%                                            trial
%                 trialType               -- trial type index (Optional)
%                 fs                      -- sampling frequency of ECoG
%                                            data
%                 T (1 x 1)               -- number of timesteps in
%                                            segment
%                 y (yDim x T)            -- neural data
%                 state (1 x T)           -- HMM state at each time
%                                            point (from the Viterbi path)
%                 p (nStates x T)       	-- defined RECURSIVELY at 
%                                            time t (1 <= t <= T) as 
%                                            the probability of the 
%                                            most probable sequence of 
%                                            length t-1 for each factor
%                                            analyzer hidden state
%                                            at time t, given the
%                                            observations from 1 to t
%                 ycs (yDim x T)          -- neural data prediction
%                 vcs (yDim x T)          -- neural data prediction
%                                            variance
%
% OPTIONAL ARGUMENTS:
%
% reconstruct - if true, no channel is 'left out' in estimating the
%               hidden Markov state Viterbi path in prediction
%               (default: false)
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  reconstruct                       = false;
  
  assignopts(who, varargin);

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

    temp                            = [];
    if (reconstruct)
      % In computing the Viterbi path, use all channels
      [temp, ~, ~]                  =...
        exactInferenceWithLL_hmm(struct('y',Yn, 'T',T),...
                                 struct('nStates',nStates,...
                                        'pi',pi, 'trans',trans,...
                                        'd',d, 'R',R,...
                                        'sharedCov', sharedCov),...
                                 'getSeq', true,...
                                 varargin{:});
    end % if (reconstruct)
    parfor i=1:yDim
      % In computing the Viterbi path, eliminate the contribution
      % of channel i
      mi                            = [1:(i-1) (i+1):yDim];
      if (reconstruct)
        state                     	= temp.state;
      else % if (~reconstruct)
        [temp2, ~, ~]               =...
          exactInferenceWithLL_hmm(struct('y',Yn(mi,:), 'T',T),...
                                   struct('nStates',nStates,...
                                          'pi',pi, 'trans',trans,...
                                          'd',d(mi,:), 'R',R(mi,mi,:),...
                                          'sharedCov', sharedCov),...
                                   'getSeq', true,...
                                   varargin{:});
        state                       = temp2.state;
      end
      
      % Taking advantage of block diagonal matrix structure
      invRmimi                    	=...
       nan(yDim-1,yDim-1,max(1,nStates*~sharedCov));
      for j=1:max(1,nStates*~sharedCov)
        if strcmp(covType, 'diagonal')
         invRmimi(:,:,j)          	= diag(diag(1./R(mi,mi,j)));
        elseif strcmp(covType, 'full')
         invRmimi(:,:,j)          	= inv(R(mi,mi,j));
        end
      end % for j=1:max(1,nStates*~sharedCov)
      temp3                         = mat2cell(invRmimi(:,:,min(state,end)),...
                                               yDim-1, yDim-1, ones(1,T));
                                    % (yDim-1)*T x (yDim-1)*T
      invBmimi                      = blkdiag(temp3{:});

      if strcmp(covType, 'full')
        temp3                       = mat2cell(R(i,mi,min(state,end)),...
                                               1, yDim-1, ones(1,T));
                                    % T x (yDim-1)*T
        Bimi                        = blkdiag(temp3{:});

        temp3                       = mat2cell(R(mi,i,min(state,end)),...
                                               yDim-1, 1, ones(1,T));
                                    % (yDim-1)*T x T
        Bmii                        = blkdiag(temp3{:});
      end % if strcmp(covType, 'full')

                                    % T x T
      Bii                           =...
       diag(ndarrayflatten(R(i,i,min(state,end))));

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