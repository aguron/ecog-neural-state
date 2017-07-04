function hmfaEngine(seqTrain, seqTest, fname, varargin)
%
% hmfaEngine(seqTrain, seqTest, fname, ...) 
%
% Extract neural trajectories using HMFA.
%
%   yDim: number of electrodes
%
% INPUTS:
%
% seqTrain    - training data structure, whose nth entry (corresponding
%               to the nth experimental trial) has fields
%                 trialId (1 x 1)         -- unique trial identifier
%                 y (# electrodes x T)    -- neural data
%                 T (1 x 1)               -- number of timesteps
% seqTest     - test data structure (same format as seqTrain)
% fname       - filename of where results are saved
%
% OPTIONAL ARGUMENTS:
%
% xDim        - state dimensionality (default: 3)
% nStates     - number of HMFA states (default: 3)
% faType      - HMFA factor analyzer(s) (with tied (0) or untied (1)
%               mean, factor loading matrix, and covariance parameters)
%               (default: [1 1 1])
% binWidth    - ECoG window width in seconds (default: 0.2)
% kernSD      - Gaussian smoothing kernel width in seconds (default: 0)

% d (yDim x nStates)              - observation mean(s)
% C (yDim x xDim x nStates)       - factor loadings (s)
% R (yDim x yDim x                - observation noise covariance(s)
%    nStates (or 1))
% pi0 (1 x nStates)               - initial HMFA start probabilities
% piPrior (1 x nStates)           - pseudo counts for HMFA start
%                                   probabilities
% trans0 (nStates x nStates)      - initial HMFA transition probabilities
% transPrior (nStates x nStates)  - pseudo counts for HMFA transition
%                                   probabilities
% stateGuess                      - initial state guesses for time points
%                                 	for training data

% prediction  - specifies whether leave-channel-out prediction
%               should be carried out (default: false)
%
% @ 2015 Akinyinka Omigbodun    aomigbod@ucsd.edu
  
  xDim                                      = 3;
  nStates                                   = 3;
  faType                                    = [1 1 1];
  binWidth                                  = 0.2;  % in sec
  kernSD                                    = 0;    % in sec

  Replicates                                = 1;
  Regularize                                = 0;
  Options                                   = [];

  d                                         = [];
  C                                         = [];
  R                                         = [];
  pi0                                       = [];
  piPrior                                   = [];
  trans0                                    = [];
  transPrior                                = [];
  stateGuess                                = [];

  learning                                  = true;
  inference                                 = true;
  prediction                                = false;

  extraOpts                                 = assignopts(who, varargin);

  if (learning)
    if ~any(faType) && (nStates > 1)
      error(['At least one of the factor analyzer parameters must ',...
             'be untied if nStates > 1']);
    end

    if ~isequal(faType,[1 1 1]) &&...
       ~isequal(faType,[1 1 0])
      fprintf('Does not support faType = [%d %d %d]\n', faType);
    end

    if (kernSD)
      fprintf('Performing presmoothing with kernSD=%g and HMFA...\n',...
              kernSD);
      % ========================
      % Presmooth data over time
      % ========================

      % Training data
      seqTrainOrig                          = seqTrain;
      for n=1:length(seqTrain)
        seqTrain(n).y                       =...
          smoother(seqTrain(n).y, kernSD, binWidth);
      end % for n=1:length(seqTrain)

      % Test data
      seqTestOrig                           = seqTest;
      for n=1:length(seqTest)
        seqTest(n).y                        =...
          smoother(seqTest(n).y, kernSD, binWidth);
      end % for n=1:length(seqTest)
    end % if (kernSD)

    % For compute efficiency, train on equal-length segments of trials
    seqTrainCut                             =...
      cutTrials(seqTrain, extraOpts{:});
    if isempty(seqTrain)
      fprintf('No segments extracted for training.\n');
    elseif isempty(seqTrainCut)
      fprintf(['WARNING: no segments extracted for training.',...
               ' Defaulting to segLength=Inf.\n']);
      seqTrainCut                           =...
        cutTrials(seqTrain, 'segLength', Inf);
    end % if isempty(seqTrainCut)

    yAll                                    = [seqTrainCut.y];
    if ~isempty(yAll)
      yDim                                  = size(yAll, 1);
    else % if isempty(yAll)
      yDim                                  = size([seqTest.y], 1);
    end

    initStatus                              =...
      ~[isempty(d) isempty(C) isempty(R),...
        isempty(pi0) isempty(piPrior) isempty(trans0) isempty(transPrior)];

    if all(~initStatus) && ~isempty(seqTrainCut)
      fprintf('Initializing parameters...\n');

      fprintf('First step: GMM...\n');
      args                                  =...
        {'Replicates', Replicates,...
         'CovType', 'diagonal', 'SharedCov', ~faType(3),...
         'Regularize', Regularize, 'Options', Options};
      if ~isempty(stateGuess)
        args                                =...
          [args, 'Start', cell2mat(stateGuess)];
      end % if ~isempty(stateGuess)
      obj                                   =...
        gmdistribution.fit(yAll', nStates, args{:});

      fprintf('Second step: MFA...\n');
      startParamsMFA.nMixComp               = nStates;
      startParamsMFA.faType                 = faType;
      startParamsMFA.Pi                     = obj.PComponents;
      startParamsMFA.d                      = obj.mu';
      for j=1:nStates
        jC                                  = j*faType(2)+(1 - faType(2));
        jR                                  = j*faType(3)+(1 - faType(3));

        startParamsMFA.C(:,:,jC)            = eye(yDim,xDim);
        startParamsMFA.R(:,:,jR)            = diag(obj.Sigma(:,:,jR));
      end % for j=1:nStates

      startParamsMFA.notes.RforceDiagonal   = true;
      estParamsMFA                          =...
        em_mfa(startParamsMFA, seqTrainCut, extraOpts{:});
      [~, seqTrainMFA]                      =...
        em_mfa(estParamsMFA, seqTrain, 'emMaxIters', 0);

      mixComp                               = cell(1,numel(seqTrainMFA));
      [mixComp{:}]                          = deal(seqTrainMFA.mixComp);
      startParams                           =...
        hmmestimate2(mixComp, mixComp, nStates, nStates,...
                     'Pseudostarts', ones(1,nStates),...
                   	 'Pseudotransitions', ones(nStates));

      % ==================================
      % Initialize state model parameters
      % ==================================
      startParams.nStates                   = estParamsMFA.nMixComp;
      startParams.faType                    = estParamsMFA.faType;
      startParams                           =...
        renamefield(startParams, {'ST', 'TR'}, {'pi', 'trans'});
      startParams                           = rmfield(startParams,'E');
      startParams.piPrior                   = ones(1, startParams.nStates);
      startParams.transPrior                = ones(startParams.nStates);

      % ========================================
      % Initialize observation model parameters
      % ========================================
      startParams.d                         = estParamsMFA.d;
      startParams.C                         = estParamsMFA.C;
      startParams.R                         = estParamsMFA.R;
    elseif ~any(~initStatus)
      fprintf('Initializing parameters with input...\n');
      % ==================================
      % Initialize state model parameters
      % ==================================
      startParams.nStates                   = nStates;
      startParams.faType                    = faType;
      startParams.pi                        = pi0;
      startParams.piPrior                   = piPrior;
      startParams.trans                     = trans0;
      startParams.transPrior                = transPrior;

      % ========================================
      % Initialize observation model parameters
      % ========================================
      startParams.d                         = d;
      startParams.C                         = C;
      startParams.R                         = R;
    else
      error(['d, C, R, pi, piPrior, trans and transPrior must',...
             ' all be specified if there is no data for training.']);
    end

    % Define parameter constraints
    startParams.notes.RforceDiagonal        = true;

    if ~isempty(seqTrainCut)
      % =====================
      % Fit model parameters
      % =====================
      fprintf('\nFitting HMFA model...\n');

      [estParams, seqTrainCut, LL, iterTime]=...
        em_hmfa(startParams, seqTrainCut, extraOpts{:});

      % Extract neural trajectories for original, unsegmented trials
      % using learned parameters
      [seqTrain, ess, LLorig]               =...
        exactInferenceWithLL_hmfa(seqTrain, estParams,...
                                  'getSeq', true,...
                                  'getLL', true,...
                                  extraOpts{:});
    else % if isempty(seqTrainCut)
      estParams                             = startParams;
    end
  else % if (~learning)
    estParams                               = loadvars(fname, 'estParams');
    if (inference)
      if ~isempty(seqTrain)
        % Extract neural trajectories for original, unsegmented trials
        % using learned parameters
        [seqTrain, ess, LLorig]            	=...
          exactInferenceWithLL_hmfa(seqTrain, estParams,...
                                    'getSeq', true,...
                                    'getLL', true,...
                                    extraOpts{:});
      end % if ~isempty(seqTrain)
    end % if (inference)
  end

  if ~isempty(seqTest) % check if there are any test trials
    % =========================================
    % Leave-channel-out prediction on test data
    % =========================================
    if prediction && estParams.notes.RforceDiagonal
      seqTest                               =...
        predict_hmfa_fast(seqTest, estParams, varargin{:});
    end

    
    % ===================================================
    % Neural trajectories and loglikelihood for test data
    % ===================================================
    if (inference)
      [seqTest, ~, LLtest]                  =...
        exactInferenceWithLL_hmfa(seqTest, estParams,...
                                  'getSeq', true,...
                                  'getLL', true,...
                                  extraOpts{:});
    end % if (inference)
  end % if ~isempty(seqTest)

  % =============
  % Save results
  % =============
  if (learning)
    vars                                  	= who;
    fprintf('Saving %s...\n', fname);
    save(fname, vars{~ismember(vars, {'yAll', 'mixComp'})});
  elseif inference || prediction
    if (inference)
      fprintf('Saving inferred latent variables in %s...\n', fname);
      if ~isempty(seqTrain)
        save(fname, 'seqTrain', 'ess', 'LLorig', '-append');
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