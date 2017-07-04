function seq = modelSample(model, T, nSeq, trialId)
%MODELSAMPLE
% 

  hmmModel.type                   = 'gauss';
  % hmmModel.modelType            = 'hmm';
  % hmmModel.emission.cpdType     = 'condGauss';
  switch(model.type)
    case 'gmm'
      hmmModel.pi                 = model.obj.PComponents;
      hmmModel.A                  = repmat(model.obj.PComponents,...
                                           [model.obj.NComponents 1]);
      nMixComp                    = model.obj.NComponents;
      % hmmModel.nstates         	= nMixComp;

      yDim                      	= size(model.obj.mu, 2);
      hmmModel.emission.d        	= yDim;
      % hmmModel.d               	= hmmModel.emission.d;

      hmmModel.emission.mu        = model.obj.mu;
      hmmModel.emission.Sigma     = zeros([yDim yDim nMixComp]);
      switch(model.obj.CovarianceType)
        case 'diagonal'
          if model.obj.SharedCovariance
           hmmModel.emission.Sigma=...
            bsxfun(@plus, hmmModel.emission.Sigma, diag(model.obj.Sigma));
          else % if ~model.obj.SharedCovariance
            for j=1:nMixComp
              hmmModel.emission.Sigma(:,:,j)...
                                  = diag(model.obj.Sigma(:,:,j));
            end % for j=1:nMixComp
          end % if model.obj.SharedCovariance
        case 'full'
          if model.obj.SharedCovariance
           hmmModel.emission.Sigma=...
            bsxfun(@plus, hmmModel.emission.Sigma, model.obj.Sigma);
          else % if ~model.obj.SharedCovariance
            for j=1:nMixComp
              hmmModel.emission.Sigma(:,:,j)...
                                  = model.obj.Sigma(:,:,j);
            end % for j=1:nMixComp
          end % if model.obj.SharedCovariance
        otherwise
          disp('Invalid covariance type');
      end % switch(model.obj.CovarianceType)
      % hmmModel.emission.nstates	= nMixComp;
      
      [y, state]                  = hmmSample(hmmModel, T, nSeq);
    case 'mfa'
      hmmModel.pi                 = model.params.Pi;
      hmmModel.A                  = repmat(model.params.Pi,...
                                           [model.params.nMixComp 1]);
      nMixComp                    = model.params.nMixComp;
      % hmmModel.nstates         	= nMixComp;

      yDim                      	= size(model.params.d, 1);
      hmmModel.emission.d        	= yDim;
      % hmmModel.d               	= hmmModel.emission.d;

      hmmModel.emission.mu        = model.params.d;
      hmmModel.emission.Sigma     = zeros([yDim yDim nMixComp]);
      faType                      = model.params.faType;
      C                           = model.params.C;
      R                           = model.params.R;
      for j=1:nMixComp
        jC                        = j*faType(2) + (1 - faType(2));
        jR                        = j*faType(3) + (1 - faType(3));
        hmmModel.emission.Sigma(:,:,j)...
                                  = C(:,:,jC)*C(:,:,jC)' + R(:,:,jR);
      end % for j=1:nMixComp
      % hmmModel.emission.nstates	= nMixComp;

      [y, state]                  = hmmSample(hmmModel, T, nSeq);
    case 'hmfa'
      hmmModel.pi                 = model.params.pi;
      hmmModel.A                  = model.params.trans;
      nStates                     = model.params.nStates;
      % hmmModel.nstates          = nStates;

      yDim                      	= size(model.params.d, 1);
      hmmModel.emission.d        	= yDim;
      % hmmModel.d               	= hmmModel.emission.d;

      hmmModel.emission.mu        = model.params.d;
      hmmModel.emission.Sigma     = zeros([yDim yDim nStates]);
      faType                      = model.params.faType;
      C                           = model.params.C;
      R                           = model.params.R;
      for j=1:nStates
        jC                        = j*faType(2) + (1 - faType(2));
        jR                        = j*faType(3) + (1 - faType(3));
        hmmModel.emission.Sigma(:,:,j)...
                                  = C(:,:,jC)*C(:,:,jC)' + R(:,:,jR);
      end % for j=1:nStates
      % hmmModel.emission.nstates	= nStates;

      [y, state]                  = hmmSample(hmmModel, T, nSeq);
    otherwise
      error(['Invalid specification of neural ',...
             'trajectory extraction method']);
  end % switch(model.type)

  if (nargin == 3)
    trialId                       = cell(1, nSeq);
    for i=1:nSeq
      trialId{i}                  =...
        sprintf(sprintf('%%0%dd', ceil(log10(nSeq+1))), i);
    end % for i=1:nSeq
  end % if (nargin == 3)
  
  if isscalar(T)
    T                             = T*ones(1,nSeq);
  end % if isscalar(T)
  
  seq                             =...
    struct('trialId', trialId,...
           'fs', num2cell(model.fs*ones(1,nSeq)),...
           'T', num2cell(T),...
           'y', y',...
           'state', state');

  switch(model.type)
    case 'gmm'
      seq                        	= renamefield(seq, 'state', 'mixComp');
    case 'mfa'
      seq                       	= renamefield(seq, 'state', 'mixComp');
      [~, temp]                   =...
        em_mfa(model.params, seq, 'emMaxIters', 0);
      seq                       	= addfield(seq, 'x', {temp.x});
    case 'hmfa'
      temp                        =...
        exactInferenceWithLL_hmfa(seq, model.params, 'getSeq', true);
      seq                       	= addfield(seq, 'x', {temp.x});
    otherwise
      error(['Invalid specification of neural ',...
             'trajectory extraction method']);
  end % switch(model.type)
end