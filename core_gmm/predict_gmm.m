function seq = predict_gmm(seq, obj, varargin)
%
% seq = predict_gmm(seq, obj, ...)
%
% Performs leave-one-channel-out prediction for GMM
%
%   yDim: number of electrodes or channels
%
% INPUTS:
%
% seq         - data structure (same format as seqTrain in GMMENGINE)
% obj         - GMM model (same format as output of MATLAB function
%               GMDISTRIBUTION.FIT)
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
%                 mixComp (1 x T)        	-- most probable mixture
%                                            component at each time point
%                 p (nMixComp x T)       	-- probability of each mixture
%                                            component at each time point
%                 ycs (yDim x T)          -- neural data prediction
%                 vcs (yDim x T)          -- neural data prediction
%                                            variance
%
% OPTIONAL ARGUMENTS:
%
% reconstruct - if true, no channel is 'left out' in estimating the
%               most probable mixture component path in prediction
%               (default: false)
%
% @ 2016 Akinyinka Omigbodun    aomigbod@ucsd.edu

  reconstruct         	= false;
  
  assignopts(who, varargin);

  yDim                  = size(obj.mu, 2);
  nMixComp              = obj.NComponents;

  Tall                  = [seq.T];
  for n=1:numel(seq)
    T                   = Tall(n);
    ycs                 = nan(yDim, T);
    vcs                 = nan(yDim, T);

    Yn                	= seq(n).y;

    parfor i=1:yDim
      % In computing the most probable component for each time point, 
      % eliminate the contribution of channel i
      mi                = [1:(i-1) (i+1):yDim];
      if (reconstruct)
        obj2           	= obj;
        mixComp        	= cluster(obj, Yn');
      else % if (~reconstruct)
        if strcmp(obj.CovType, 'diagonal')
          obj2         	= gmdistribution(obj.mu(:,mi),...
                                         obj.Sigma(1,mi,:),...
                                         obj.PComponents);
        elseif strcmp(obj.CovType, 'full')
          obj2         	= gmdistribution(obj.mu(:,mi),...
                                         obj.Sigma(mi,mi,:),...
                                         obj.PComponents);
        end
        mixComp        	= cluster(obj2, Yn(mi,:)');
      end

      % Taking advantage of block diagonal matrix structure
      invSigma        	=...
       nan(yDim-1,yDim-1,max(1,nMixComp*~obj.SharedCov));
      for j=1:max(1,nMixComp*~obj.SharedCov)
        if strcmp(obj.CovType, 'diagonal')
         invSigma(:,:,j)= diag(1./obj2.Sigma(:,:,j));
        elseif strcmp(obj.CovType, 'full')
         invSigma(:,:,j)= inv(obj2.Sigma(:,:,j));
        end
      end % for j=1:max(1,nMixComp*~obj.SharedCov)
      temp            	= mat2cell(invSigma(:,:,min(mixComp,end)),...
                                   yDim-1, yDim-1, ones(1,T));
                        % (yDim-1)*T x (yDim-1)*T
      invBmimi         	= blkdiag(temp{:});


      if strcmp(obj.CovType, 'full')
        temp          	= mat2cell(obj.Sigma(i,mi,min(mixComp,end)),...
                                   1, yDim-1, ones(1,T));
                        % T x (yDim-1)*T
        Bimi          	= blkdiag(temp{:});
      end

      if strcmp(obj.CovType, 'full')
        temp            = mat2cell(obj.Sigma(mi,i,min(mixComp,end)),...
                                   yDim-1, 1, ones(1,T));
                        % (yDim-1)*T x T
        Bmii          	= blkdiag(temp{:});
      end

                        % T x T
      if strcmp(obj.CovType, 'diagonal')
        Bii            	= diag(ndarrayflatten(obj.Sigma(1,i,min(mixComp,end))));
      elseif strcmp(obj.CovType, 'full')
        Bii            	= diag(ndarrayflatten(obj.Sigma(i,i,min(mixComp,end))));
      end

      if strcmp(obj.CovType, 'diagonal')
        % 1 x T
        ycs(i,:)        = obj.mu(mixComp,i)';
      elseif strcmp(obj.CovType, 'full')
                        % yDim-1 x T
        difmi          	= Yn(mi,:) - obj.mu(mixComp,mi)';

        % 1 x T
        ycs(i,:)        =...
          obj.mu(mixComp,i)' + (Bimi*invBmimi*difmi(:))';
      end

      if strcmp(obj.CovType, 'diagonal')
        % 1 x T
        vcs(i,:)        = diag(Bii)';
      elseif strcmp(obj.CovType, 'full')
        % 1 x T
        vcs(i,:)        = diag(Bii - Bimi*invBmimi*Bmii)';
      end
    end % parfor i=1:yDim
    seq(n).ycs          = ycs;
    seq(n).vcs          = vcs;
    fprintf('Cross-validation complete for trial %d\r', n);
  end % for n=1:numel(seq)
  fprintf('\n');
end