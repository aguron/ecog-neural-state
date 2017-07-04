function result = neuralstate(runIdx, dat, varargin)
%
% result = neuralstate(runIdx, dat, ...)
%
% Prepares data and calls functions for extracting neural trajectories.
%
% INPUTS:
%
% runIdx        - results files will be saved in mat_results/runXXX, where
%                 XXX is runIdx
% dat           - structure whose nth entry (corresponding to the nth
%                 experimental trial) has fields
%                   trialId       -- unique trial identifier
%                   trialType     -- index (OPTIONAL)
%                   fs            -- sampling frequency of ECoG data
%
%                   ECoG          -- matrix of voltage activity across all
%                                    electrodes. Each row corresponds to an
%                                    electrode. Each column corresponds to
%                                    a (1/fs) sec timestep.
%
%                  	OR
%
%                   y (yDim x T)  -- neural data
%                   T (1 x 1)     -- number of timesteps
%
% OUTPUTS:
%
% result        - structure containing all variables saved in
%                 mat_results/runXXX/ if 'numFolds' is 0.
%                 Else, the structure is empty.
%               
% OPTIONAL ARGUMENTS:
%
% prediction    - specifies whether leave-channel-out prediction
%                 should be carried out (default: false)
% predOverride  - overwrite current leave-channel-out prediction
%                 (default: false)
% method        - method for extracting neural trajectories
%                 'mhmfa', 'mhmm', 'chmfa', 'phmfa', 'hmfa', 'mfa', 'gmm',
%                 'gpfa' (default), 'fa', 'ppca', 'pca'
% binWidth      - ECoG window width in sec (default: 0.2)
% numFolds      - number of cross-validation folds (default: 0)
%                 0 indicates no cross-validation, i.e. train
%                 on all trials only
% fracTrainData - fraction of training data used for each cross-validation
%                 fold
% xDim          - state dimensionality (default: 8)
% nStates       - number of HMFA states (default: 3)
% nMixComp      - number of MFA/GMM/MHMFA/MHMM mixture components
%                 (default: 3)
% faType        - CHMFA/PHMFA/HMFA/MFA factor analyzer(s) (with tied (0)
%                 or untied (1) mean, factor loading matrix, and covariance
%                 parameters) (default: [1 1 1])
% CovType       - GMM covariance type: 'full' or 'diagonal' (default)
% mixCompGuess 	- initial component guesses for time points or MHMFA/MHMM
%                 mixture components
% stateGuess  	- initial state guesses for time points
% partition     - ECoG channel assignment to coupled/parallel
%                 Markov processes
%
% @ 2009 Byron Yu         byronyu@stanford.edu
%        John Cunningham  jcunnin@stanford.edu

  prediction                        = false;
  predOverride                      = false;
  inferOverride                     = false;
  method                            = 'gpfa';
  binWidth                          = 0.2; % in sec
  numFolds                          = 0;
  fracTrainData                     = 1;
  xDim                              = 8;
  nStates                           = 3;
  nMixComp                          = 3;
  faType                            = [1 1 1];
  CovType                           = 'diagonal';
  SharedCov                         = false;
  mixCompGuess                      = [];
  stateGuess                       	= [];
  partition                         = [];
  nAttempts                         = Inf;
  extraOpts                         = assignopts(who, varargin);

  if ismember(method, {'chmfa','phmfa'})
    if isempty(partition) || (numel(partition) ~= size(dat(1).ECoG,1))
      fprintf(['Error: Invalid partition of channels into ',...
               'Markov processes.  Exiting.\n']);
      result                        = [];
      return
    end

    if ~isequal(1:max(partition),unique(partition))
      error(['All the integers between 1 and max(partition)',...
             ' must be present in partition']);
    end % if ~isequal(1:max(partition),unique(partition))
    
    nProcesses                      = max(partition); % Markov Processes
    if (nProcesses < 2)
      fprintf(['Error: There must be at least 2',...
               'Markov processes.  Exiting.\n']);
      result                        = [];
      return
    end % if (nProcesses < 2)
  else % if ~ismember(method, {'chmfa','phmfa'})
    if (numel(xDim) ~= 1) ||...
       (numel(nStates) ~= 1) ||...
       (iscell(faType) && (numel(faType) ~= 1))
      error('Invalid parameter input specification');
    end
  end

  fprintf('\n---------------------------------------\n');
  if ~isdir('mat_results')
    mkdir('mat_results');
  end % if ~isdir('mat_results')
  % Make a directory for this runIdx if it doesn't already exist
  runDir                            =...
    sprintf('mat_results/run%03d/binWidth_%g', runIdx, binWidth);
  if isdir(runDir)
    fprintf('Using existing directory %s...\n', runDir);
  else
    fprintf('Making directory %s...\n', runDir);
    mkdir(runDir);
  end % if isdir(runDir)

  % Obtain RMS power
  seq                               = getSeqRMS(dat, binWidth);
  if isempty(seq)
    fprintf('Error: No valid trials.  Exiting.\n');
    result                          = [];
    return
  end
  % Set cross-validation folds
  N                                 = numel(seq);
  % Randomly reorder trials before partitioning into
  % training and test sets
  rng('default')
  tr                                = randperm(N);
  if (numFolds)
    if isfield(seq, 'trialType')
      cvfs                         	=...
        preparecvfs([seq(tr).trialType], numFolds);
    else % if ~isfield(seq, 'trialType')
      fdiv                         	= floor(linspace(1, N+1, numFolds+1));
    end
  end % if (numFolds)
  for cvf=0:numFolds
    if cvf == 0
      fprintf('\n===== Training on all data =====\n');
    else
      fprintf('\n===== Cross-validation fold %d of %d =====\n',...
              cvf, numFolds);
    end

    % Specify filename where results will be saved
    switch(method)
      case 'mhmfa'
        faTypeSpec                  = 'tu'; % t - tied; u - untied
        fname                       =...
          sprintf('%s/%s_nMixComp%02d_xDim%02d_nStates%02d_MFV%c%c%c',...
                  runDir, method, nMixComp, xDim, nStates,...
                  faTypeSpec(faType+1));
      case 'mhmm'
        fname                       =...
          sprintf('%s/%s_nMixComp%02d_nStates%02d_%s',...
                  runDir, method, nMixComp, nStates, CovType);
      case {'chmfa', 'phmfa'}
        fname                       = sprintf('%s/%s', runDir, method);
        if (numel(xDim) ~= 1) && (numel(xDim) ~= nProcesses)
          fprintf('Invalid xDim parameter specification');
          result                   	= [];
          return
        else
          fname                     = sprintf('%s_xDim%02d',...
                                              fname, xDim(1));
          for sp=xDim(2:end)
            fname                   = sprintf('%s_%02d', fname, sp);
          end % for sp=xDim(2:end)
        end
        if (numel(nStates) ~= 1) && (numel(nStates) ~= nProcesses)
          fprintf('Invalid nStates parameter specification');
          result                   	= [];
          return
        else
          fname                     = sprintf('%s_nStates%02d',...
                                              fname, nStates(1));
          for sp=nStates(2:end)
            fname                   = sprintf('%s_%02d', fname, sp);
          end % for sp=nStates(2:end)
        end
        faTypeSpec                  = 'tu'; % t - tied; u - untied
        if iscell(faType)
          if (numel(faType) ~= nProcesses)
            fprintf('Invalid faType parameter specification');
            result                 	= [];
            return
          else
            fname                   = sprintf('%s_MFV%c%c%c', fname,...
                                              faTypeSpec(faType{1}+1));
            for sp=faType(2:end)
              fname                 = sprintf('%s_%c%c%c', fname, sp{1});
            end % for sp=faType(2:end)
          end
        else
          fname                     = sprintf('%s_MFV%c%c%c', fname,...
                                              faTypeSpec(faType+1));
        end
      case 'hmfa'
        faTypeSpec                  = 'tu'; % t - tied; u - untied
        fname                       =...
          sprintf('%s/%s_xDim%02d_nStates%02d_MFV%c%c%c',...
                  runDir, method, xDim, nStates, faTypeSpec(faType+1));
      case 'mfa'
        faTypeSpec                  = 'tu'; % t - tied; u - untied
        fname                       =...
          sprintf('%s/%s_xDim%02d_nMixComp%02d_MFV%c%c%c',...
                  runDir, method, xDim, nMixComp, faTypeSpec(faType+1));
      case {'gmm', 'hmm'}
       switch(method)
        case 'gmm'
         fname                      =...
          sprintf('%s/%s_nMixComp%02d_%s',...
                  runDir, method, nMixComp, CovType);
        case 'hmm'
         fname                     	=...
          sprintf('%s/%s_nStates%02d_%s',...
                  runDir, method, nStates, CovType);
       end % switch(method)
       if (SharedCov)
        fname                      	= sprintf('%s_tied',fname);
       end % if (SharedCov)
      case {'gpfa', 'fa', 'ppca', 'pca'}
        fname                       =...
          sprintf('%s/%s_xDim%02d', runDir, method, xDim);
      otherwise
        error(['Invalid specification of neural ',...
               'trajectory extraction method']);
    end % switch(method)
    if (cvf > 0)
      fname                         = sprintf('%s_cv%02dof%02d',...
                                              fname, cvf, numFolds);
    end % if (cvf > 0)
    learning                        = true;
    inference                       = true;
    if exist([fname '.mat'], 'file')
      if (~inferOverride)
        inference                   = false;
      end % if (~inferOverride)
     
      if exist('predFlag','var') && (~prediction)
       prediction                   = true;
      end % if exist('predFlag','var') && (~prediction)
      if ((inference) || ((prediction) && (cvf > 0))) &&...
         ismember(method, {'chmfa', 'phmfa', 'hmfa', 'mfa', 'hmm', 'gmm'})
        seqTest                     = loadvars(fname, 'seqTest');

        if (isfield(seqTest,'ycsOrth01') || isfield(seqTest,'ycs'))...
           && (~predOverride)
          fprintf(['%s.mat already exists with leave-one-channel-out ',...
                   'prediction.'], fname);
          if (prediction)
           predFlag                 = true;
          end % if (prediction)
          prediction                = false;
          if (~inferOverride)
            fprintf('  Skipping...\n');
            continue
          else
            fprintf('  Inferring latent variables...\n');
            learning                = false;
          end
        else
          learning                 	= false;
        end
      else
        fprintf('%s.mat already exists.  Skipping...\n', fname);
        continue
      end
    end % if exist([fname '.mat'], 'file')

    if (learning)
      testMask                      = false(1, N);
      if (cvf > 0)
        % Balance trial types if trialType is a field of seq
        if isfield(seq, 'trialType')
          Ntrain                    = sum(cvfs{cvf});
          testMask                  = ~cvfs{cvf};
        else % if ~isfield(seq, 'trialType')
          Ntrain                   	= N - numel(fdiv(cvf):fdiv(cvf+1)-1);
          testMask(fdiv(cvf):fdiv(cvf+1)-1)...
                                    = true;
        end
      else % if (cvf == 0)
        Ntrain                     	= N;
      end
      trainMask                    	= ~testMask;
      
      % test set
      testTrialIdx                 	= tr(testMask);
      seqTest                      	= seq(testTrialIdx);
      
      % specified fraction of training set
      trainTrialIdx                	= tr(trainMask);
      Ntrain                        = ceil(Ntrain*fracTrainData);
      if isfield(seq, 'trialType')
        trainTrialIdx               =...
        trainTrialIdx(balanceselectedtrials([seq(trainTrialIdx).trialType],...
                                            Ntrain));
      else % if ~isfield(seq, 'trialType')
        trainTrialIdx              	= trainTrialIdx(1:Ntrain);
      end
      seqTrain                     	= seq(trainTrialIdx);
    else % if (~learning)
      seqTrain                      = loadvars(fname, 'seqTrain');
      if ismember(method, {'chmfa','phmfa'})
        partition                  	= loadvars(fname, 'partition');
      end % if ismember(method, {'chmfa','phmfa'})
    end

    % Check if training data covariance is full rank
    yAll                            = [seqTrain.y];

    if ~isempty(yAll)
      if ismember(method, {'chmfa','phmfa'})
        for p=1:nProcesses
          yDim                        = size(yAll(partition==p,:), 1);
          if rank(cov(yAll(partition==p,:)')) < yDim
            fprintf(['ERROR: Observation covariance matrix for ',...
                     'process %d is rank deficient.\n'], p);
            fprintf(['Possible causes: repeated units, not enough ',...
                     'observations.\n']);
            fprintf('Exiting...\n');
            return
          end
        end % for p=1:nProcesses
      else
        yDim                          = size(yAll, 1);
        if rank(cov(yAll')) < yDim
          fprintf('ERROR: Observation covariance matrix is rank deficient.\n');
          fprintf('Possible causes: repeated units, not enough observations.\n');
          fprintf('Exiting...\n');
          return
        end
      end
    end % if ~isempty(yAll)
    
    fprintf('Number of training trials: %d\n', length(seqTrain));
    fprintf('Number of test trials: %d\n', length(seqTest));

    % If doing cross-validation, don't use private noise variance floor.  
    if (cvf == 1)
      extraOpts                     = [extraOpts, 'minVarFrac', -Inf];
    end % if (cvf == 1)

    % The following does the heavy lifting.
    done                            = false;
    attempt                         = 0;
    while(~done)
     try
      switch(method)
        case 'mhmfa'
          fprintf('Latent space dimensionality: %d\n', xDim);
          fprintf('Number of hidden Markov states: %d\n', nStates);
          fprintf('Number of MHMFA mixture components: %d\n', nMixComp);
          if ~isempty(stateGuess)
            if ~eq(numel(stateGuess),numel(seq))
              error('Invalid specification of stateGuess');
            end % if ~eq(numel(stateGuess),numel(seq))
            extraOpts               =...
              [extraOpts, 'stateGuess', stateGuess(trainTrialIdx)];
          end % if ~isempty(stateGuess)
          if ~isempty(mixCompGuess)
            if ~eq(numel(mixCompGuess),numel(seq))
              error('Invalid specification of mixCompGuess');
            end % if ~eq(numel(mixCompGuess),numel(seq))
            extraOpts               =...
              [extraOpts, 'mixCompGuess', mixCompGuess(trainTrialIdx)];
          end % if ~isempty(mixCompGuess)
          mhmfaEngine(seqTrain, seqTest, fname,...
                     'xDim', xDim, 'nStates', nStates, 'faType', faType,...
                     'nMixComp', nMixComp, 'binWidth', binWidth,...
                     'learning', learning, 'inference', inference,...
                     'prediction', prediction,...
                     extraOpts{:});
        case 'mhmm'
          fprintf('Number of hidden Markov states: %d\n', nStates);
          fprintf('Number of MHMM mixture components: %d\n', nMixComp);
          if ~isempty(stateGuess)
            if ~eq(numel(stateGuess),numel(seq))
              error('Invalid specification of stateGuess');
            end % if ~eq(numel(stateGuess),numel(seq))
            extraOpts               =...
              [extraOpts, 'stateGuess', stateGuess(trainTrialIdx)];
          end % if ~isempty(stateGuess)
          if ~isempty(mixCompGuess)
            if ~eq(numel(mixCompGuess),numel(seq))
              error('Invalid specification of mixCompGuess');
            end % if ~eq(numel(mixCompGuess),numel(seq))
            extraOpts               =...
              [extraOpts, 'mixCompGuess', mixCompGuess(trainTrialIdx)];
          end % if ~isempty(mixCompGuess)
          extraOpts                 =...
            [extraOpts, 'xDim', xDim, 'faType', faType];
          mhmmEngine(seqTrain, seqTest, fname,...
                     'nStates', nStates, 'nMixComp', nMixComp,...
                     'CovType', CovType, 'binWidth', binWidth,...
                     'learning', learning, 'inference', inference,...
                     'prediction', prediction,...
                     extraOpts{:});
        case 'chmfa'
          for p=1:nProcesses
            fprintf('Process %d\n', p);
            fprintf('Latent space dimensionality: %d\n', xDim(min(end,p)));
            fprintf('Number of hidden Markov states: %d\n',...
                    nStates(min(end,p)));
          end % for p=1:nProcesses
          if ~isempty(stateGuess)
            if ~eq(numel(stateGuess),numel(seq))
              error('Invalid specification of stateGuess');
            end % if ~eq(numel(stateGuess),numel(seq))
            extraOpts               =...
              [extraOpts, 'stateGuess', stateGuess(trainTrialIdx)];
          end % if ~isempty(stateGuess)
          chmfaEngine(seqTrain, seqTest, partition, fname,...
                      'xDim', xDim, 'nStates', nStates, 'faType', faType,...
                      'binWidth', binWidth, 'learning', learning,...
                      'inference', inference, 'prediction', prediction,...
                      extraOpts{:});
        case 'phmfa'
          for p=1:nProcesses
            fprintf('Process %d\n', p);
            fprintf('Latent space dimensionality: %d\n', xDim(min(end,p)));
            fprintf('Number of hidden Markov states: %d\n',...
                    nStates(min(end,p)));
          end % for p=1:nProcesses
          if ~isempty(stateGuess)
            if ~eq(numel(stateGuess),numel(seq))
              error('Invalid specification of stateGuess');
            end % if ~eq(numel(stateGuess),numel(seq))
            extraOpts               =...
              [extraOpts, 'stateGuess', stateGuess(trainTrialIdx)];
          end % if ~isempty(stateGuess)
          phmfaEngine(seqTrain, seqTest, partition, fname,...
                      'xDim', xDim, 'nStates', nStates, 'faType', faType,...
                      'binWidth', binWidth, 'prediction', prediction,...
                     'predOnly', predOnly, extraOpts{:});
        case 'hmfa'
          fprintf('Latent space dimensionality: %d\n', xDim);
          fprintf('Number of hidden Markov states: %d\n', nStates);
          if ~isempty(stateGuess)
            if ~eq(numel(stateGuess),numel(seq))
              error('Invalid specification of stateGuess');
            end % if ~eq(numel(stateGuess),numel(seq))
            extraOpts               =...
              [extraOpts, 'stateGuess', stateGuess(trainTrialIdx)];
          end % if ~isempty(stateGuess)
          hmfaEngine(seqTrain, seqTest, fname,...
                     'xDim', xDim, 'nStates', nStates, 'faType', faType,...
                     'binWidth', binWidth, 'learning', learning,...
                     'inference', inference, 'prediction', prediction,...
                     extraOpts{:});
        case 'mfa'
          fprintf('Latent space dimensionality: %d\n', xDim);
          fprintf('Number of factor analyzer mixture components: %d\n',...
                  nMixComp);
          if ~isempty(mixCompGuess)
            extraOpts               =...
              [extraOpts, 'mixCompGuess', mixCompGuess(trainTrialIdx)];
          end % if ~isempty(mixCompGuess)
          mfaEngine(seqTrain, seqTest, fname,...
                    'xDim', xDim, 'nMixComp', nMixComp, 'faType', faType,...
                    'binWidth', binWidth, 'learning', learning,...
                     'inference', inference, 'prediction', prediction,...
                     extraOpts{:});
        case 'hmm'
          fprintf('Number of hidden Markov states: %d\n', nStates);
          if ~isempty(stateGuess)
            if ~eq(numel(stateGuess),numel(seq))
              error('Invalid specification of stateGuess');
            end % if ~eq(numel(stateGuess),numel(seq))
            extraOpts               =...
              [extraOpts, 'stateGuess', stateGuess(trainTrialIdx)];
          end % if ~isempty(stateGuess)
          hmmEngine(seqTrain, seqTest, fname,...
                    'nStates', nStates, 'CovType', CovType,... 
                    'binWidth', binWidth, 'learning', learning,...
                    'inference', inference, 'prediction', prediction,...
                    'SharedCov', SharedCov, extraOpts{:});
        case 'gmm'
          fprintf('Number of mixture components: %d\n', nMixComp);
          if ~isempty(mixCompGuess)
            extraOpts               =...
              [extraOpts, 'mixCompGuess', mixCompGuess(trainTrialIdx)];
          end % if ~isempty(mixCompGuess)
          gmmEngine(seqTrain, seqTest, fname,...
                    'nMixComp', nMixComp, 'CovType', CovType,...
                    'binWidth', binWidth, 'learning', learning,...
                     'inference', inference, 'prediction', prediction,...
                     'SharedCov', SharedCov, extraOpts{:});
        case 'gpfa'
          fprintf('Latent space dimensionality: %d\n', xDim);
          gpfaEngine2(seqTrain, seqTest, fname,...
                      'xDim', xDim, 'binWidth', binWidth, extraOpts{:});
        case {'fa', 'ppca', 'pca'}
          fprintf('Latent space dimensionality: %d\n', xDim);
          twoStageEngine2(seqTrain, seqTest, fname,...
                          'typ', method, 'xDim', xDim,...
                          'binWidth', binWidth, extraOpts{:});
        otherwise
          error(['Invalid specification of neural ',...
                 'trajectory extraction method']);
      end % switch(method)

      if exist([fname '.mat'], 'file')
        save(fname, 'method', 'cvf', 'extraOpts', '-append');
      end % if exist([fname '.mat'], 'file')
      
      done                          = true;
     catch err
      attempt                       = attempt + 1;
      displayerror(err);
      fprintf('To discontinue model fit, set done = true\n');
      programcontrol
      if (done) || (attempt == nAttempts)
        rethrow(err)
      end % if (done) || (attempt == nAttempts)
     end % try
    end % while(~done)
  end % for cvf=0:numFolds

  result                            = [];
  if (nargout == 1)
    if (cvf > 0)
      fname                         = rstrtok(fname,'_');
    end % if (cvf > 0)
    if exist([fname '.mat'], 'file')
      result                      	= load(fname);
      if ismember(method, {'chmfa', 'phmfa', 'hmfa', 'mfa'})
        if ~isfield(result, 'AIC')
         [result.AIC result.nParams]= aic(method,...
                                          result.estParams,...
                                          result.LLorig);
         AIC                       	= result.AIC;
         nParams                   	= result.nParams;
         save(fname, 'AIC', 'nParams', '-append');
        end % if ~isfield(result, 'AIC')

        if ~isfield(result, 'BIC')
          result.BIC               	= bic(method,...
                                          result.estParams,...
                                          result.LLorig,...
                                          sum([result.seqTrain.T]));
          BIC                     	= result.BIC;
          save(fname, 'BIC', '-append');
        end % if ~isfield(result, 'BIC')
      end % if ismember(method, {'chmfa', 'phmfa', 'hmfa', 'mfa'})
    end % if exist([fname '.mat'], 'file')
  end % if (nargout == 1)
end
