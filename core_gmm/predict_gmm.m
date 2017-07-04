function seq = predict_gmm(seq, obj, varargin)
%
% seq = predict_gmm(seq, obj, ...)
%
% Performs leave-electrode-out prediction for GMM.
%
% INPUTS:
%
% seq       - test data structure
% obj       - GMM model fit to training data
%
% OUTPUTS:
%
% seq      	- test data structure with new fields
%               ycs (yDim x T)
%               vcs (yDim x T)
%
% @ 2016 Akinyinka Omigbodun    aomigbod@ucsd.edu

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
      % eliminate the contribution of electrode i
      mi                = [1:(i-1) (i+1):yDim];

      if strcmp(obj.CovType, 'diagonal')
        obj2            = gmdistribution(obj.mu(:,mi),...
                                         obj.Sigma(1,mi,:),...
                                         obj.PComponents);
      elseif strcmp(obj.CovType, 'full')
       	obj2            = gmdistribution(obj.mu(:,mi),...
                                         obj.Sigma(mi,mi,:),...
                                         obj.PComponents);
      end
      mixComp           = cluster(obj2, Yn(mi,:)');

      % Taking advantage of block diagonal matrix structure
      invSigma        	= nan(yDim-1,yDim-1,nMixComp);
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
        Bii            	= diag(matflat(obj.Sigma(1,i,min(mixComp,end))));
      elseif strcmp(obj.CovType, 'full')
        Bii            	= diag(matflat(obj.Sigma(i,i,min(mixComp,end))));
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