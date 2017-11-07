function mhmfaEngine(seqTrain, seqTest, fname, varargin)
%
% mhmfaEngine(seqTrain, seqTest, fname, ...) 
%
% Model fitting and inference with MHMFA
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
% xDim       	- state dimensionality (default: 3)
% nMixComp   	- number of MHMFA mixture components (default: 3)
% nStates    	- number of Markov (HMFA) states (default: 3)
% faType     	- component HMFA factor analyzers specification
%               (with tied (0) or untied (1) mean, factor loading matrix,
%               and covariance parameters) (default: [1 1 1])
%
% binWidth  	- ECoG window width in seconds (default: 0.2)
% kernSD      - Gaussian smoothing kernel width in seconds. 0 corresponds
%               to no smoothing (default: 0)
%
% d {nMixComp x 1} (yDim x nStates)               - cell array of component
%                                                   HMFA observation means
% C {nMixComp x 1} (yDim x xDim x nStates (or 1))	- cell array of component
%                                                   HMFA factor loadings
% R {nMixComp x 1} (yDim x yDim x nStates (or 1))	- cell array of component
%                                                   HMFA observation noise
%                                                   covariances
%
% pi0 {nMixComp x 1} (1 x nStates)               	- cell array of component
%                                                   HMFA initial start
%                                                   probabilities
% piPrior {nMixComp x 1} (1 x nStates)           	- cell array of pseudo
%                                                   counts for component
%                                                   HMFA start
%                                                   probabilities
% trans0 {nMixComp x 1} (nStates x nStates)     	- cell array of component
%                                                   HMFA initial transition
%                                                   probabilities
% transPrior {nMixComp x 1} (nStates x nStates) 	- cell array of pseudo
%                                                   counts for component
%                                                   HMFA transition
%                                                   probabilities
% Pi (1 x nMixComp)                               - component HMFA priors
%
% stateGuess                                      - cell array of initial
%                                                   state guesses for time
%                                                 	points of trials for
%                                                  	training data
%
% mixCompGuess                                    - vector of initial
%                                                   component HMFA guesses
%                                                   for training trials.
%                                                   Use an index of 0 if
%                                                   there is no guess for a
%                                                   trial; an index between
%                                                   1 and nMixComp
%                                                   (inclusive), if the
%                                                   trial is to be used in
%                                                   initialization; and an
%                                                   index between -nMixComp
%                                                   and -1 (inclusive), if
%                                                  	the trial is not to be
%                                                   used in model fitting
%	removeOutliers                                 	- if true, outlier trials
%                                                   identified during the
%                                                   initialization stage
%                                                   are eliminated from
%                                                   model fitting
%                                                   (default: true)
%
% learning                                        - indicates whether model
%                                                   fitting should be
%                                                   carried out
%                                                   (default: true)
% inference                                       - indicates whether
%                                                   latent variables should
%                                                   be inferred for seqTest
%                                                   (default: true)
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu
  
  xDim                                      = 3;
  nMixComp                                  = 3;
  nStates                                   = 3;
  faType                                    = [1 1 1];

  binWidth                                  = 0.2;
  kernSD                                    = 0;
  
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

  extraOpts                                 = assignopts(who, varargin);
  
  if (learning)
    if ~any(faType) && (nStates > 1)
      error(['At least one of the factor analyzer parameters must ',...
             'be untied if nStates > 1']);
    end

    if ~isequal(faType,[1 1 1]) &&...
       ~isequal(faType,[1 1 0]) &&...
       ~isequal(faType,[1 0 0])
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
      fprintf('No segments for training.\n');
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
      startParams(1).faType                	= faType;
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
    % ===================================================
    % Latent variables and loglikelihood for test data
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
  elseif (inference)
    fprintf('Saving inferred latent variables in %s...\n', fname);
    if ~isempty(seqTrain)
      save(fname, 'seqTrain', 'ess', 'LLorig', '-append');
    end % if ~isempty(seqTrain)
    if ~isempty(seqTest)
      save(fname, 'seqTest', 'LLtest', '-append');
    end % if ~isempty(seqTest)
    save(fname, 'inference', '-append');
  end
end