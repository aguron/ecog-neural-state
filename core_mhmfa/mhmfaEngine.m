function mhmfaEngine(seqTrain, seqTest, fname, varargin)
%
% mhmfaEngine(seqTrain, seqTest, fname, ...) 
%
% Extract neural trajectories using MHMFA.
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
% nMixComp    - number of MHMFA mixture components (default: 3)
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
% Pi (1 x nMixComp)               - MHMFA component priors

% stateGuess                      - initial state guesses for time points
%                                 	for training data
% mixCompGuess                    - initial MHMFA component guesses
%                                   for training data: use an index of
%                                   0 if there is no guess for a trial
%                                   mixComp; an index between 1 and
%                                   nMixComp inclusive, if the mixComp
%                                   guess is to be used in the
%                                   initialization; and an index between
%                                   -nMixComp and -1 inclusive, if the
%                                   mixComp guess is not to be used in the
%                                   initialization

% learning    - specifies whether model fitting should be carried out
%               (default: true)
% inference   - specifies whether the latent variables should be
%               inferred for seqTest (default: true)
% prediction  - specifies whether leave-channel-out prediction
%               should be carried out on seqTest (default: false)
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu
  
  xDim                                      = 3;
  nStates                                   = 3;
  faType                                    = [1 1 1];
  nMixComp                                  = 3;
  binWidth                                  = 0.2;  % in sec
  kernSD                                    = 0;    % in sec
  
  d                                         = [];
  C                                         = [];
  R                                         = [];
  pi0                                       = [];
  piPrior                                   = [];
  trans0                                    = [];
  transPrior                                = [];
  Pi                                        = [];

  stateGuess                                = [];
  mixCompGuess                              = [];
  removeOutliers                            = true;

  learning                                  = true;
  inference                                 = true;
  prediction                                = false;

  extraOpts                                 = assignopts(who, varargin);
  
  if (learning)
    if ~any(faType) && (nStates > 1)
      error(['At least one of the factor analyzer parameters must ',...
             'be untied if nStates > 1']);
    end % if ~any(faType) && (nStates > 1)

    if ~isequal(faType,[1 1 1]) &&...
       ~isequal(faType,[1 1 0])
      fprintf('Does not support faType = [%d %d %d]\n', faType);
    end

    if (kernSD)
      fprintf('Performing presmoothing with kernSD=%g and MHMFA...\n',...
              kernSD);
      % ========================
      % Presmooth data over time
      % ========================

      % Training data
      seqTrainOrig                          = seqTrain;
      for n=1:numel(seqTrain)
        seqTrain(n).y                       =...
          smoother(seqTrain(n).y, kernSD, binWidth);
      end % for n=1:numel(seqTrain)

      % Test data
      seqTestOrig                           = seqTest;
      for n=1:numel(seqTest)
        seqTest(n).y                        =...
          smoother(seqTest(n).y, kernSD, binWidth);
      end % for n=1:numel(seqTest)
    end % if (kernSD)

    % Resegment trials for training
    [seqTrainCut, resegTrlGuess]           	=...
      resegmenttrials(seqTrain,...
                      'method', 'mhmfa',...
                      'stateGuess', stateGuess,...
                      'mixCompGuess', mixCompGuess,...
                      extraOpts{:});
    if isempty(seqTrain)
      fprintf('No segments extracted for training.\n');
    elseif isempty(seqTrainCut)
      fprintf(['WARNING: no segments extracted for training.',...
               ' Defaulting to segLength=Inf.\n']);
      [seqTrainCut, resegTrlGuess]         	=...
        resegmenttrials(seqTrain,...
                        'method', 'mhmfa',...
                        'stateGuess', stateGuess,...
                        'mixCompGuess', mixCompGuess,...
                        'segLength', Inf);
    end

    yAll                                    = [seqTrainCut.y];
    if ~isempty(yAll)
      yDim                                  = size(yAll, 1);
    else % if isempty(yAll)
      yDim                                  = size([seqTest.y], 1);
    end

    initStatus                              =...
      ~[isempty(d) isempty(C) isempty(R),...
        isempty(pi0) isempty(piPrior),...
        isempty(trans0) isempty(transPrior),...
        isempty(Pi)];

    if all(~initStatus) && ~isempty(seqTrainCut)
      if removeOutliers
        argsOut                             = cell(3,1);
      else % if ~removeOutliers
        argsOut                             = cell(2,1);
      end
      [argsOut{:}]                          =...
        init_mhmfa(seqTrainCut, fname,...
                   'xDim', xDim, 'nStates', nStates,...
                   'faType', faType, 'nMixComp', nMixComp,...
                   'stateGuess', resegTrlGuess.state,...
                   'mixCompGuess', resegTrlGuess.mixComp,...
                   extraOpts{:});
      startParams                           = argsOut{1};
      resegTrlGuess.mixComp               	= argsOut{2};
      if removeOutliers
        outliers                            = argsOut{3};
      end % if removeOutliers

      if (numunique(resegTrlGuess.mixComp(abs(resegTrlGuess.mixComp)>0))...
          ~= nMixComp)
        error('Invalid mixCompGuess specification');
      end
    elseif ~any(~initStatus)
      fprintf('Initializing parameters with input...\n');
      % ==================================
      % Initialize state model parameters
      % ==================================
      startParams(1).nStates               	= nStates;
      startParams(1).faType                	= faType;
      startParams(1).nMixComp              	= nMixComp;

      startParams(nMixComp).Pi              = num2cell(Pi);
      [startParams.Pi]                      =...
        deal(startParams(nMixComp).Pi{:});
      startParams(nMixComp).pi              = class2cell(pi0);
      [startParams.pi]                      =...
        deal(startParams(nMixComp).pi{:});
      startParams(nMixComp).piPrior        	= class2cell(piPrior);
      [startParams.piPrior]                	=...
        deal(startParams(nMixComp).piPrior{:});
      startParams(nMixComp).trans           = class2cell(trans0);
      [startParams.trans]                   =...
        deal(startParams(nMixComp).trans{:});
      startParams(nMixComp).transPrior     	= class2cell(transPrior);
      [startParams.transPrior]             	=...
        deal(startParams(nMixComp).transPrior{:});

      % ========================================
      % Initialize observation model parameters
      % ========================================
      startParams(nMixComp).d               = class2cell(d);
      [startParams.d]                       =...
        deal(startParams(nMixComp).d{:});
      startParams(nMixComp).C               = class2cell(C);
      [startParams.C]                       =...
        deal(startParams(nMixComp).C{:});
      startParams(nMixComp).R               = class2cell(R);
      [startParams.R]                       =...
        deal(startParams(nMixComp).R{:});

      % Define parameter constraints
      startParams(1).notes.RforceDiagonal  	= true;
    else
      error(['d, C, R, Pi, pi, piPrior, trans and transPrior must',...
             ' all be specified if there is no data for training.']);
    end

    if ~isempty(seqTrainCut)
      % =====================
      % Fit model parameters
      % =====================
      fprintf('\nFitting MHMFA model...\n');

      [estParams, seqTrainCut, LL, iterTime]=...
        em_mhmfa(startParams, seqTrainCut,...
                 'outliers', outliers, extraOpts{:});

      % Inference for original, unsegmented trials using learned parameters
      [seqTrain, ess, LLorig]               =...
        exactInferenceWithLL_mhmfa(seqTrain, estParams,...
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
        % Inference for original, unsegmented trials using learned
        % parameters
        [seqTrain, ess, LLorig]            	=...
          exactInferenceWithLL_mhmfa(seqTrain, estParams,...
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
    if prediction && estParams(1).notes.RforceDiagonal
      
    end

    % ===================================================
    % Neural trajectories and loglikelihood for test data
    % ===================================================
    if (inference)
      [seqTest, ~, LLtest]                  =...
        exactInferenceWithLL_mhmfa(seqTest, estParams,...
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
    save(fname, vars{~ismember(vars, {'yAll'})});
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