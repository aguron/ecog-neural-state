function result = neuralstate(runIdx, dat, varargin)
%
% result = neuralstate(runIdx, dat, ...)
%
% Prepares data and calls functions for model fitting, inference,
%	and/or prediction
%
% INPUTS:
%
% runIdx        - results files will be saved in
%                 mat_results/runXXX/binWidth_O...O/, where
%                 XXX is runIdx and O...O is binWidth
% dat           - structure whose n-th entry (corresponding to the n-th
%                 experimental trial) has fields
%                   trialId       -- unique trial identifier
%                   trialType     -- trial type index (Optional)
%                   fs            -- sampling frequency of ECoG data
%            EITHER:
%                   ECoG          -- matrix of voltage activity across all
%                                    electrodes. Each row corresponds to an
%                                    electrode. Each column corresponds to
%                                    a (1/fs) sec timestep.
%                OR:
%                   T (1 x 1)     -- number of timesteps
%                   y (yDim x T)  -- neural data
%
% OUTPUTS:
%
% result        - structure containing all variables saved in a particular
%                 results file in mat_results/runXXX/binWidth_O...O/
%                 if 'numFolds' is 0. Else, the structure is empty.
%               
% OPTIONAL ARGUMENTS:
%
% prediction    - specifies whether leave-channel-out prediction
%                 should be carried out (default: false)
% predOverride  - overwrite current leave-channel-out prediction
%                 (default: false)
% inferOverride - overwrite current inferred latent variables
%                 (default: false)
% method        - method for model fitting, inference, and/or prediction:
%                 'mhmfa', 'mhmm', 'hmfa' (default), 'mfa', 'hmm', 'gmm'
% binWidth      - ECoG window width in seconds (please also see RUNIDX
%                 above) (default: 0.2)
% numFolds      - number of cross-validation folds (default: 0)
%                 0 indicates no cross-validation, i.e. train
%                 on all trials only
% fracTrainData - fraction of training data used for each cross-validation
%                 fold or of all the data when there is no cross-validation
%                 (default: 1)
% xDim          - state dimensionality (default: 3)
% nStates       - number of Markov states (default: 3)
% nMixComp      - number of MHMFA/MHMM/MFA/GMM mixture components
%                 (default: 3)
% faType        - MHMFA/HMFA/MFA factor analyzers specification
%                 (with tied (0) or untied (1) mean, factor loading matrix,
%                 and covariance parameters) (default: [1 1 1])
% CovType       - MHMM/HMM/GMM covariance type: 'full' or
%                 'diagonal' (default)
% SharedCov     - MHMM (mixture components)/HMM/GMM covariance tied (true)
%                 or untied (false) (default: false)
% mixCompGuess 	- initial mixture component guesses for trials of
%                 MHMFA/MHMM (vector) or time points of trials for
%                 MFA/GMM (cell array) 
% stateGuess  	- cell array of initial state guesses for time points
%                 of trials for MHMFA/MHMM/HMFA/HMM
% nAttempts     - number of model fitting attempts
%
% Code adapted from neuralTraj.m by Byron Yu and John Cunningham.
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  prediction                        = false;
  predOverride                      = false;
  inferOverride                     = false;
  method                            = 'hmfa';
  binWidth                          = 0.2;
  numFolds                          = 0;
  fracTrainData                     = 1;
  xDim                              = 3;
  nStates                           = 3;
  nMixComp                          = 3;
  faType                            = [1 1 1];
  CovType                           = 'diagonal';
  SharedCov                         = false;
  mixCompGuess                      = [];
  stateGuess                       	= [];
  nAttempts                         = Inf;
  extraOpts                         = assignopts(who, varargin);

  fprintf('\n---------------------------------------\n');
  if ~isdir('mat_results')
    mkdir('mat_results');
  end % if ~isdir('mat_results')
  % Make a directory for this runIdx if it doesn't already exist
  runDir                            =...
    sprintf('mat_results/run%03d/binWidth_%g', runIdx, binWidth);
  if isdir(runDir)
    fprintf('Using existing directory %s...\n', runDir);
  else % if ~isdir(runDir)
    fprintf('Making directory %s...\n', runDir);
    mkdir(runDir);
  end

  result                            = [];

  % Obtain RMS power
  seq                               = getSeqRMS(dat, binWidth);
  if isempty(seq)
    fprintf('Error: No valid trials.  Exiting.\n');
    result                          = [];
    return
  end % if isempty(seq)
  % Set cross-validation folds
  N                                 = numel(seq);
  % Randomly reorder trials before partitioning into training and test sets
  rng('default')
  tr                                = randperm(N);
  if (numFolds)
    % The proportion of trial types (if specified) is similar across
    % cross-validation folds
    if isfield(seq, 'trialType')
      cvfs                         	= preparecvfs([seq(tr).trialType],...
                                                  numFolds);
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
        if (SharedCov)
         fname                    	= sprintf('%s_tied',fname);
        end % if (SharedCov)
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
      otherwise
        error(['Invalid specification of neural',... 
               'state inference method']);
    end % switch(method)
    if (cvf > 0)
      fname                         = sprintf('%s_cv%02dof%02d',...
                                              fname, cvf, numFolds);
    end % if (cvf > 0)
    learning                        = true;
    inference                       = true;
    if exist('predFlag','var') && (~prediction)
      prediction                    = true;
    end % if exist('predFlag','var') && (~prediction)
    if exist([fname '.mat'], 'file')
      if (~inferOverride)
        inference                   = false;
      end % if (~inferOverride)

      if ((inference) || ((prediction) && (cvf > 0))) &&...
         ismember(method, {'hmfa', 'mfa', 'hmm', 'gmm'})
        seqTest                     = loadvars(fname, 'seqTest');

        if isfield(seqTest,'ycs') && (~predOverride)
          fprintf(['%s.mat already exists with leave-one-channel-out ',...
                   'prediction.'], fname);
          if (prediction)
           predFlag                 = true;
          end % if (prediction)
          prediction                = false;
          if (~inferOverride)
            fprintf('  Skipping...\n');
            continue
          else % if (inferOverride)
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
          % Vector entry is true in a cross-validation fold when
          % it belongs to the training set, and false otherwise
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
        % Attempt to equalize the numbers of trial types (if specified)
        % for the selected trials
        trainTrialIdx               =...
        trainTrialIdx(balanceselectedtrials([seq(trainTrialIdx).trialType],...
                                            Ntrain));
      else % if ~isfield(seq, 'trialType')
        trainTrialIdx              	= trainTrialIdx(1:Ntrain);
      end
      seqTrain                     	= seq(trainTrialIdx);
    else % if (~learning)
      seqTrain                      = loadvars(fname, 'seqTrain');
    end

    % Check if training data covariance is full rank
    yAll                            = [seqTrain.y];

    if ~isempty(yAll)
      yDim                          = size(yAll, 1);
      if rank(cov(yAll')) < yDim
        fprintf('ERROR: Observation covariance matrix is rank deficient.\n');
        fprintf('Possible causes: repeated units, not enough observations.\n');
        fprintf('Exiting...\n');
        return
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
          fprintf('Number of MHMFA mixture components: %d\n\n', nMixComp);
          fprintf('Number of hidden Markov states: %d\n', nStates);
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
                     extraOpts{:});
        case 'mhmm'
          fprintf('Number of MHMM mixture components: %d\n\n', nMixComp);
          fprintf('Number of hidden Markov states: %d\n', nStates);
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
                     'SharedCov', SharedCov,...
                     extraOpts{:});
        case 'hmfa'
          fprintf('Latent space dimensionality: %d\n', xDim);
          fprintf('Number of hidden Markov states: %d\n\n', nStates);
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
          fprintf('Number of factor analyzer mixture components: %d\n\n',...
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
          fprintf('Number of hidden Markov states: %d\n\n', nStates);
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
          fprintf('Number of mixture components: %d\n\n', nMixComp);
          if ~isempty(mixCompGuess)
            extraOpts               =...
              [extraOpts, 'mixCompGuess', mixCompGuess(trainTrialIdx)];
          end % if ~isempty(mixCompGuess)
          gmmEngine(seqTrain, seqTest, fname,...
                    'nMixComp', nMixComp, 'CovType', CovType,...
                    'binWidth', binWidth, 'learning', learning,...
                     'inference', inference, 'prediction', prediction,...
                     'SharedCov', SharedCov, extraOpts{:});
        otherwise
          error('Invalid specification of neural state method');
      end % switch(method)

      if exist([fname '.mat'], 'file')
        save(fname, 'method', 'cvf', 'attempt', 'extraOpts', '-append');
      end % if exist([fname '.mat'], 'file')
      
      done                          = true;
     catch err
      attempt                       = attempt + 1;
      displayerror(err);
      fprintf('Attempt %d unsuccessful\n\n', attempt);
      % Please see PROGRAMCONTROL in misc-functionality/debugging
      programcontrol
      % To discontinue model fitting, inference, and/or prediction:
      %   set done = true
      if (done) || (attempt == nAttempts)
        rethrow(err)
      end % if (done) || (attempt == nAttempts)
     end % try
    end % while(~done)

    if (nargout == 1) && exist([fname '.mat'], 'file')
      result(cvf+1).model           = load(fname);
    end
  end % for cvf=0:numFolds
end