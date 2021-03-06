function mfaEngine(seqTrain, seqTest, fname, varargin)
%
% mfaEngine(seqTrain, seqTest, fname, ...) 
%
% Model fitting and inference with MFA
%
%   yDim: number of electrodes or channels
%
% INPUTS:
%
% seqTrain    - training data structure, whose n-th entry (corresponding
%               to the n-th experimental trial) has fields
%                 trialId           -- unique trial identifier
%                 trialType (1 x 1)	-- trial type index (Optional)
%                 fs (1 x 1)       	-- sampling frequency of ECoG data
%                 T (1 x 1)         -- number of timesteps
%                 y (yDim x T)      -- neural data
% seqTest     - test data structure (same format as seqTrain)
% fname       - model filename
%
% OPTIONAL ARGUMENTS:
%
% xDim                              - state dimensionality (default: 3)
% nMixComp                          - number of MFA mixture components
%                                     (default: 3)
% faType                            - factor analyzers specification
%                                     (with tied (0) or untied (1) mean,
%                                     factor loading matrix, and covariance
%                                     parameters) (default: [1 1 1])
%
% binWidth                          - ECoG window width in seconds
%                                     (default: 0.2)
% kernSD                            - Gaussian smoothing kernel width in
%                                    	seconds. 0 corresponds to no
%                                     smoothing (default: 0)
%
% Replicates                       ]  parameters for
% Regularize                       ]- MATLAB function
% Options                          ]  GMDISTRIBUTION.FIT
%
% d (yDim x nMixComp)             	- MFA observation means
% C (yDim x xDim x nMixComp (or 1))	- MFA factor loadings
% R (yDim x yDim x nMixComp (or 1))	- MFA observation noise covariances
%
% Pi (1 x nMixComp)                	- component factor analyzer priors
%
% mixCompGuess                    	- cell array of initial component
%                                     guesses for time points of trials for
%                                     training data
%
% learning                          - indicates whether model fitting
%                                     should be carried out (default: true)
% inference                       	- indicates whether latent variables
%                                     should be inferred for seqTest
%                                     (default: true)
% prediction                        - specifies whether
%                                     leave-one-channel-out prediction
%                                    	should be carried out
%                                     (default: false)
%
% @ 2016 Akinyinka Omigbodun    aomigbod@ucsd.edu

  xDim                              = 3;
  nMixComp                          = 3;
  faType                            = [1 1 1];

  binWidth                          = 0.2;
  kernSD                            = 0;
  
  Replicates                        = 1;
  Regularize                        = 0;
  Options                           = [];

  d                                 = [];
  C                                 = [];
  R                                 = [];

  Pi                                = [];
  
  mixCompGuess                      = [];

  learning                          = true;
  inference                         = true;
  prediction                        = false;

  extraOpts                         = assignopts(who, varargin);

  if (learning)
    if ~any(faType) && (nMixComp > 1)
      error(['At least one of the factor analyzer parameters must ',...
             'be untied if nMixComp > 1']);
    end % if ~any(faType) && (nMixComp > 1)

    if ~isequal(faType,[1 1 1]) &&...
       ~isequal(faType,[1 1 0]) &&...
       ~isequal(faType,[1 0 0])
     fprintf('Does not support faType = [%d %d %d]\n',faType);
    end

    if (kernSD)
      fprintf('Performing presmoothing with kernSD=%g and MFA...\n',...
              kernSD);
      % ========================
      % Presmooth data over time
      % ========================

      % Training data
      seqTrainOrig                  = seqTrain;
      for n=1:length(seqTrain)
        seqTrain(n).y               =...
         smoother(seqTrain(n).y, kernSD, binWidth);
      end % for n=1:length(seqTrain)

      % Test data
      seqTestOrig                   = seqTest;
      for n=1:length(seqTest)
        seqTest(n).y                =...
         smoother(seqTest(n).y, kernSD, binWidth);
      end % for n=1:length(seqTest)
    end % if (kernSD)

    % Resegment trials for training
    [seqTrainCut, resegTrlGuess]    =...
      resegmenttrials(seqTrain,...
                      'method', 'mfa',...
                      'mixCompGuess', mixCompGuess,...
                      extraOpts{:});
    if isempty(seqTrain)
      fprintf('No segments for training.\n');
    elseif isempty(seqTrainCut)
      fprintf(['WARNING: no segments extracted for training.',...
               ' Defaulting to segLength=Inf.\n']);
      [seqTrainCut, resegTrlGuess]	=...
        resegmenttrials(seqTrain,...
                        'method', 'mfa',...
                        'mixCompGuess', mixCompGuess,...
                        'segLength', Inf);
    end

    yAll                            = [seqTrainCut.y];
    if ~isempty(yAll)
      yDim                          = size(yAll, 1);
    else % if isempty(yAll)
      yDim                          = size([seqTest.y], 1);
    end

    initStatus                      =...
      ~[isempty(d) isempty(C) isempty(R) isempty(Pi)];
    if all(~initStatus) && ~isempty(seqTrainCut)
      fprintf('Initializing parameters using GMM...\n');

      args                          =...
        {'Replicates', Replicates,...
         'CovType', 'diagonal', 'SharedCov', ~faType(3),...
         'Regularize', Regularize, 'Options', Options};
      if ~isempty(resegTrlGuess.mixComp)
        args                        =...
         [args, 'Start', cell2mat(resegTrlGuess.mixComp)];
      end % if ~isempty(resegTrlGuess.mixComp)
      obj                           =...
       gmdistribution.fit(yAll',nMixComp,args{:});
    end % if all(~initStatus) && ~isempty(seqTrainCut)

    % ==================================
    % Initialize mixture model parameters
    % ==================================
    startParams.nMixComp            = nMixComp;
    if exist('obj','var')
      startParams.Pi                = obj.PComponents;
    elseif ~any(~initStatus)
      fprintf('Initializing parameters with input...\n');
      startParams.Pi                = Pi;
    else
      error(['d, C, R, and Pi must all be specified',...
             ' if there is no data for training.']);
    end

    % ========================================
    % Initialize observation model parameters
    % ========================================
    if exist('obj','var')
      startParams.d                 = obj.mu';
      for j=1:nMixComp
        jC                          = j*faType(2) + (1 - faType(2));
        jR                          = j*faType(3) + (1 - faType(3));

        startParams.C(:,:,jC)       = eye(yDim,xDim);
        startParams.R(:,:,jR)       = diag(obj.Sigma(:,:,jR));
      end % for j=1:nMixComp
    elseif ~any(~initStatus)
      startParams.d                 = d;
      startParams.C                 = C;
      startParams.R                 = R;
    end

    % Define parameter constraints
    startParams.notes.RforceDiagonal= true;
   	startParams.faType              = faType;                      

    if ~isempty(seqTrainCut)
      % =====================
      % Fit model parameters
      % =====================
      fprintf('\nFitting MFA model...\n');

      [estParams, seqTrainCut, LL, iterTime]...
                                   	=...
        em_mfa(startParams, seqTrainCut, extraOpts{:});

      % Extract neural trajectories for original, unsegmented trials
      % using learned parameters
      [~, seqTrain, LLorig]         =...
        em_mfa(estParams, seqTrain, 'emMaxIters', 0);
    else % if isempty(seqTrainCut)
      estParams                     = startParams;
    end
  else % if (~learning)
    estParams                       = loadvars(fname, 'estParams');
    if (inference)
      if ~isempty(seqTrain)
        % Extract neural trajectories for original, unsegmented trials
        % using learned parameters
        [~, seqTrain, LLorig]       =...
          em_mfa(estParams, seqTrain, 'emMaxIters', 0);
      end % if ~isempty(seqTrain)
    end % if (inference)
  end
  
  if ~isempty(seqTest) % check if there are any test trials
    % =============================================
    % Leave-one-channel-out prediction on test data
    % =============================================
    if prediction && estParams.notes.RforceDiagonal
      seqTest                       =...
        predict_mfa_fast(seqTest, estParams, varargin{:});
    end

    % ===================================================
    % Neural trajectories and loglikelihood for test data
    % ===================================================
    if (inference)
      [~, seqTest, LLtest]          =...
        em_mfa(estParams, seqTest, 'emMaxIters', 0);
    end % if (inference)
  end % if ~isempty(seqTest)

  % =============
  % Save results
  % =============
  if (learning)
    vars                            = who;
    fprintf('Saving %s...\n', fname);
    save(fname, vars{~ismember(vars, {'yAll'})});
  elseif inference || prediction
    if (inference)
      fprintf('Saving inferred latent variables in %s...\n', fname);
      if ~isempty(seqTrain)
        save(fname, 'seqTrain', 'LLorig', '-append');
      end % if ~isempty(seqTrain)
      if ~isempty(seqTest)
        save(fname, 'seqTest', 'LLtest', '-append');
      end % if ~isempty(seqTest)
      save(fname, 'inference', '-append');
    end % if (inference)

    if (prediction) && ~isempty(seqTest)
      fprintf('Saving leave-one-channel-out prediction in %s...\n', fname);
      if (inference)
        save(fname, 'prediction', '-append');
      else % if (~inference)
        save(fname, 'seqTest', 'prediction', '-append');
      end
    end
  end
end